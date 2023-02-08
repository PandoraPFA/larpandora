/**
 *  @file   larpandora/LArPandoraEventBuilding/LArPandoraTrackCreation_module.cc
 *
 *  @brief  module for lar pandora track creation
 */

#include "larpandoracontent/LArObjects/LArPfoObjects.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"

#include <memory>

namespace lar_pandora {

  class LArPandoraTrackCreation : public art::EDProducer {
  public:
    explicit LArPandoraTrackCreation(fhicl::ParameterSet const& pset);

    LArPandoraTrackCreation(LArPandoraTrackCreation const&) = delete;
    LArPandoraTrackCreation(LArPandoraTrackCreation&&) = delete;
    LArPandoraTrackCreation& operator=(LArPandoraTrackCreation const&) = delete;
    LArPandoraTrackCreation& operator=(LArPandoraTrackCreation&&) = delete;

    void produce(art::Event& evt) override;

  private:
    /**
     *  @brief Build a recob::Track object
     *
     *  @param id the id code for the track
     *  @param trackStateVector the vector of trajectory points for this track
     */
    recob::Track BuildTrack(const int id,
                            const lar_content::LArTrackStateVector& trackStateVector) const;

    std::string m_pfParticleLabel;       ///< The pf particle label
    unsigned int m_minTrajectoryPoints;  ///< The minimum number of trajectory points
    unsigned int m_slidingFitHalfWindow; ///< The sliding fit half window
    bool m_useAllParticles;              ///< Build a recob::Track for every recob::PFParticle
  };

  DEFINE_ART_MODULE(LArPandoraTrackCreation)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"

#include "art/Persistency/Common/PtrMaker.h"

#include "canvas/Utilities/InputTag.h"

#include "larcore/Geometry/Geometry.h"

#include "lardata/Utilities/AssociationUtil.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandora/LArPandoraInterface/Detectors/GetDetectorType.h"
#include "larpandora/LArPandoraInterface/Detectors/LArPandoraDetectorType.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <iostream>

namespace lar_pandora {

  LArPandoraTrackCreation::LArPandoraTrackCreation(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    , m_pfParticleLabel(pset.get<std::string>("PFParticleLabel"))
    , m_minTrajectoryPoints(pset.get<unsigned int>("MinTrajectoryPoints", 2))
    , m_slidingFitHalfWindow(pset.get<unsigned int>("SlidingFitHalfWindow", 20))
    , m_useAllParticles(pset.get<bool>("UseAllParticles", false))
  {
    produces<std::vector<recob::Track>>();
    produces<art::Assns<recob::PFParticle, recob::Track>>();
    produces<art::Assns<recob::Track, recob::Hit>>();
    produces<art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>>();

    if (m_minTrajectoryPoints < 2)
      throw cet::exception("LArPandoraTrackCreation")
        << "MinTrajectoryPoints should not be smaller than 2!";
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraTrackCreation::produce(art::Event& evt)
  {
    std::unique_ptr<std::vector<recob::Track>> outputTracks(new std::vector<recob::Track>);
    std::unique_ptr<art::Assns<recob::PFParticle, recob::Track>> outputParticlesToTracks(
      new art::Assns<recob::PFParticle, recob::Track>);
    std::unique_ptr<art::Assns<recob::Track, recob::Hit>> outputTracksToHits(
      new art::Assns<recob::Track, recob::Hit>);
    std::unique_ptr<art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>>
      outputTracksToHitsWithMeta(new art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>);

    std::unique_ptr<LArPandoraDetectorType> detType(detector_functions::GetDetectorType());
    // 'wirePitchW` is here used only to provide length scale for binning hits and performing sliding/local linear fits.
    const float wirePitchW(detType->WirePitchW());

    int trackCounter(0);
    const art::PtrMaker<recob::Track> makeTrackPtr(evt);

    // Organise inputs
    PFParticleVector pfParticleVector, extraPfParticleVector;
    PFParticlesToSpacePoints pfParticlesToSpacePoints;
    PFParticlesToClusters pfParticlesToClusters;
    LArPandoraHelper::CollectPFParticles(
      evt, m_pfParticleLabel, pfParticleVector, pfParticlesToSpacePoints);
    LArPandoraHelper::CollectPFParticles(
      evt, m_pfParticleLabel, extraPfParticleVector, pfParticlesToClusters);

    VertexVector vertexVector;
    PFParticlesToVertices pfParticlesToVertices;
    LArPandoraHelper::CollectVertices(evt, m_pfParticleLabel, vertexVector, pfParticlesToVertices);

    for (const art::Ptr<recob::PFParticle> pPFParticle : pfParticleVector) {
      // Select track-like pfparticles
      if (!m_useAllParticles && !LArPandoraHelper::IsTrack(pPFParticle)) continue;

      // Obtain associated spacepoints
      PFParticlesToSpacePoints::const_iterator particleToSpacePointIter(
        pfParticlesToSpacePoints.find(pPFParticle));

      if (pfParticlesToSpacePoints.end() == particleToSpacePointIter) {
        mf::LogDebug("LArPandoraTrackCreation") << "No spacepoints associated to particle ";
        continue;
      }

      // Obtain associated clusters
      PFParticlesToClusters::const_iterator particleToClustersIter(
        pfParticlesToClusters.find(pPFParticle));

      if (pfParticlesToClusters.end() == particleToClustersIter) {
        mf::LogDebug("LArPandoraShowerCreation") << "No clusters associated to particle ";
        continue;
      }

      // Obtain associated vertex
      PFParticlesToVertices::const_iterator particleToVertexIter(
        pfParticlesToVertices.find(pPFParticle));

      if ((pfParticlesToVertices.end() == particleToVertexIter) ||
          (1 != particleToVertexIter->second.size())) {
        mf::LogDebug("LArPandoraTrackCreation") << "Unexpected number of vertices for particle ";
        continue;
      }

      // Copy information into expected pandora form
      pandora::CartesianPointVector cartesianPointVector;
      for (const art::Ptr<recob::SpacePoint> spacePoint : particleToSpacePointIter->second)
        cartesianPointVector.emplace_back(pandora::CartesianVector(
          spacePoint->XYZ()[0], spacePoint->XYZ()[1], spacePoint->XYZ()[2]));

      double vertexXYZ[3] = {0., 0., 0.};
      particleToVertexIter->second.front()->XYZ(vertexXYZ);
      const pandora::CartesianVector vertexPosition(vertexXYZ[0], vertexXYZ[1], vertexXYZ[2]);

      // Call pandora "fast" track fitter
      lar_content::LArTrackStateVector trackStateVector;
      pandora::IntVector indexVector;
      try {
        lar_content::LArPfoHelper::GetSlidingFitTrajectory(cartesianPointVector,
                                                           vertexPosition,
                                                           m_slidingFitHalfWindow,
                                                           wirePitchW,
                                                           trackStateVector,
                                                           &indexVector);
      }
      catch (const pandora::StatusCodeException&) {
        mf::LogDebug("LArPandoraTrackCreation") << "Unable to extract sliding fit trajectory";
        continue;
      }

      if (trackStateVector.size() < m_minTrajectoryPoints) {
        mf::LogDebug("LArPandoraTrackCreation")
          << "Insufficient input trajectory points to build track: " << trackStateVector.size();
        continue;
      }

      HitVector hitsFromSpacePoints, hitsFromClusters, hitsInParticle;
      HitSet hitsInParticleSet;

      LArPandoraHelper::GetAssociatedHits(evt,
                                          m_pfParticleLabel,
                                          particleToSpacePointIter->second,
                                          hitsFromSpacePoints,
                                          &indexVector);
      LArPandoraHelper::GetAssociatedHits(
        evt, m_pfParticleLabel, particleToClustersIter->second, hitsFromClusters);
      //ATTN: hits ordered from space points if available, rest added at the end
      for (unsigned int hitIndex = 0; hitIndex < hitsFromSpacePoints.size(); hitIndex++) {
        hitsInParticle.push_back(hitsFromSpacePoints.at(hitIndex));
        (void)hitsInParticleSet.insert(hitsFromSpacePoints.at(hitIndex));
      }

      for (unsigned int hitIndex = 0; hitIndex < hitsFromClusters.size(); hitIndex++) {
        if (hitsInParticleSet.count(hitsFromClusters.at(hitIndex)) == 0)
          hitsInParticle.push_back(hitsFromClusters.at(hitIndex));
      }

      // Add invalid points at the end of the vector, so that the number of the trajectory points is the same as the number of hits
      if (trackStateVector.size() > hitsFromSpacePoints.size()) {
        throw cet::exception("LArPandoraTrackCreation")
          << "trackStateVector.size() is greater than hitsFromSpacePoints.size()";
      }
      const unsigned int nInvalidPoints = hitsInParticle.size() - trackStateVector.size();
      for (unsigned int i = 0; i < nInvalidPoints; ++i) {
        trackStateVector.push_back(lar_content::LArTrackState(
          pandora::CartesianVector(util::kBogusF, util::kBogusF, util::kBogusF),
          pandora::CartesianVector(util::kBogusF, util::kBogusF, util::kBogusF),
          nullptr));
      }

      // Output objects
      outputTracks->emplace_back(
        LArPandoraTrackCreation::BuildTrack(trackCounter++, trackStateVector));
      art::Ptr<recob::Track> pTrack(makeTrackPtr(outputTracks->size() - 1));

      // Output associations, after output objects are in place
      util::CreateAssn(*this, evt, pTrack, pPFParticle, *(outputParticlesToTracks.get()));
      util::CreateAssn(
        *this, evt, *(outputTracks.get()), hitsInParticle, *(outputTracksToHits.get()));

      //ATTN: metadata added with index from space points if available, null for others
      for (unsigned int hitIndex = 0; hitIndex < hitsInParticle.size(); hitIndex++) {
        const art::Ptr<recob::Hit> pHit(hitsInParticle.at(hitIndex));
        const int index((hitIndex < hitsFromSpacePoints.size()) ? hitIndex :
                                                                  std::numeric_limits<int>::max());
        recob::TrackHitMeta metadata(index, -std::numeric_limits<double>::max());
        outputTracksToHitsWithMeta->addSingle(pTrack, pHit, metadata);
      }
    }

    mf::LogDebug("LArPandoraTrackCreation")
      << "Number of new tracks: " << outputTracks->size() << std::endl;

    evt.put(std::move(outputTracks));
    evt.put(std::move(outputTracksToHits));
    evt.put(std::move(outputTracksToHitsWithMeta));
    evt.put(std::move(outputParticlesToTracks));
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  recob::Track LArPandoraTrackCreation::BuildTrack(
    const int id,
    const lar_content::LArTrackStateVector& trackStateVector) const
  {
    if (trackStateVector.empty())
      throw cet::exception("LArPandoraTrackCreation")
        << "BuildTrack - No input trajectory points provided ";

    recob::tracking::Positions_t xyz;
    recob::tracking::Momenta_t pxpypz;
    recob::TrackTrajectory::Flags_t flags;

    for (const lar_content::LArTrackState& trackState : trackStateVector) {
      xyz.emplace_back(recob::tracking::Point_t(trackState.GetPosition().GetX(),
                                                trackState.GetPosition().GetY(),
                                                trackState.GetPosition().GetZ()));
      pxpypz.emplace_back(recob::tracking::Vector_t(trackState.GetDirection().GetX(),
                                                    trackState.GetDirection().GetY(),
                                                    trackState.GetDirection().GetZ()));
      // Set flag NoPoint if point has bogus coordinates, otherwise use clean flag set
      if (std::fabs(trackState.GetPosition().GetX() - util::kBogusF) <
            std::numeric_limits<float>::epsilon() &&
          std::fabs(trackState.GetPosition().GetY() - util::kBogusF) <
            std::numeric_limits<float>::epsilon() &&
          std::fabs(trackState.GetPosition().GetZ() - util::kBogusF) <
            std::numeric_limits<float>::epsilon()) {
        flags.emplace_back(recob::TrajectoryPointFlags(recob::TrajectoryPointFlags::InvalidHitIndex,
                                                       recob::TrajectoryPointFlagTraits::NoPoint));
      }
      else {
        flags.emplace_back(recob::TrajectoryPointFlags());
      }
    }

    // note from gc: eventually we should produce a TrackTrajectory, not a Track with empty covariance matrix and bogus chi2, etc.
    return recob::Track(
      recob::TrackTrajectory(std::move(xyz), std::move(pxpypz), std::move(flags), false),
      util::kBogusI,
      util::kBogusF,
      util::kBogusI,
      recob::tracking::SMatrixSym55(),
      recob::tracking::SMatrixSym55(),
      id);
  }

} // namespace lar_pandora
