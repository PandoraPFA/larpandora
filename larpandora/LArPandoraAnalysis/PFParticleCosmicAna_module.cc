/**
 *  @file   larpandora/LArPandoraAnalysis/PFParticleCosmicAna_module.cc
 *
 *  @brief  Analysis module for created particles
 */

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "TTree.h"

#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora {

  /**
 *  @brief  PFParticleCosmicAna class
 */
  class PFParticleCosmicAna : public art::EDAnalyzer {
  public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    PFParticleCosmicAna(fhicl::ParameterSet const& pset);

    void beginJob();
    void endJob();
    void analyze(const art::Event& evt);
    void reconfigure(fhicl::ParameterSet const& pset);

  private:
    /**
     *  @brief Fill event-level variables using input maps between reconstructed objects
     *
     *  @param  recoParticlesToHits  mapping from particles to hits
     *  @param  recoParticlesToTracks  mapping from particles to tracks
     *  @param  recoTracksToCosmicTags  mapping from tracks to cosmic tags
     */
    void FillRecoTree(const PFParticlesToHits& recoParticlesToHits,
                      const PFParticlesToTracks& recoParticlesToTracks,
                      const TracksToCosmicTags& recoTracksToCosmicTags);

    /**
     *  @brief Fill track-level variables using input maps between reconstructed objects
     *
     *  @param  hitVector  input vector of reconstructed hits
     *  @param  trueHitsToParticles  mapping between true hits and particles
     *  @param  recoHitsToParticles  mapping between reconstructed hits and particles
     *  @param  particlesToTruth  mapping between MC particles and MC truth
     *  @param  particlesToTracks  mapping between reconstructed particles and tracks
     *  @param  tracksToCosmicTags  mapping between reconstructed tracks and cosmic tags
     */
    void FillTrueTree(const HitVector& hitVector,
                      const HitsToMCParticles& trueHitsToParticles,
                      const HitsToPFParticles& recoHitsToParticles,
                      const MCParticlesToMCTruth& particlesToTruth,
                      const PFParticlesToTracks& particlesToTracks,
                      const TracksToCosmicTags& tracksToCosmicTags);

    /**
     *  @brief Get cosmic score for a PFParticle using track-level information
     *
     *  @param  particle  input reconstructed particle
     *  @param  recoParticlesToTracks  mapping between reconstructed particles and tracks
     *  @param  recoTracksToCosmicTags  mapping between reconstructed tracks and cosmic tags
     */
    float GetCosmicScore(const art::Ptr<recob::PFParticle> particle,
                         const PFParticlesToTracks& recoParticlesToTracks,
                         const TracksToCosmicTags& recoTracksToCosmicTags) const;

    TTree* m_pRecoTree; ///<
    TTree* m_pTrueTree; ///<

    int m_run;   ///<
    int m_event; ///<
    int m_index; ///<

    int m_self;          ///<
    int m_pdgCode;       ///<
    int m_isTrackLike;   ///<
    int m_isPrimary;     ///<
    float m_cosmicScore; ///<
    int m_nTracks;       ///<
    int m_nHits;         ///<

    float m_trackVtxX;       ///<
    float m_trackVtxY;       ///<
    float m_trackVtxZ;       ///<
    float m_trackEndX;       ///<
    float m_trackEndY;       ///<
    float m_trackEndZ;       ///<
    float m_trackVtxDirX;    ///<
    float m_trackVtxDirY;    ///<
    float m_trackVtxDirZ;    ///<
    float m_trackEndDirX;    ///<
    float m_trackEndDirY;    ///<
    float m_trackEndDirZ;    ///<
    float m_trackLength;     ///<
    float m_trackWidthX;     ///<
    float m_trackWidthY;     ///<
    float m_trackWidthZ;     ///<
    float m_trackVtxDeltaYZ; ///<
    float m_trackEndDeltaYZ; ///<

    int m_trackVtxContained; ///<
    int m_trackEndContained; ///<

    int m_nNeutrinoHits;
    int m_nNeutrinoHitsFullyTagged;
    int m_nNeutrinoHitsSemiTagged;
    int m_nNeutrinoHitsNotTagged;
    int m_nNeutrinoHitsNotReconstructed;
    int m_nNeutrinoHitsReconstructed;

    int m_nCosmicHits;
    int m_nCosmicHitsFullyTagged;
    int m_nCosmicHitsSemiTagged;
    int m_nCosmicHitsNotTagged;
    int m_nCosmicHitsNotReconstructed;
    int m_nCosmicHitsReconstructed;

    std::string m_hitfinderLabel;   ///<
    std::string m_trackfitLabel;    ///<
    std::string m_particleLabel;    ///<
    std::string m_cosmicLabel;      ///<
    std::string m_geantModuleLabel; ///<

    bool m_useDaughterPFParticles; ///<
    bool m_useDaughterMCParticles; ///<

    double m_cosmicContainmentCut; ///<
  };

  DEFINE_ART_MODULE(PFParticleCosmicAna)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include <iostream>

namespace lar_pandora {

  PFParticleCosmicAna::PFParticleCosmicAna(fhicl::ParameterSet const& pset) : art::EDAnalyzer(pset)
  {
    m_cosmicLabel = pset.get<std::string>("CosmicTagModule", "cosmictagger");
    m_particleLabel = pset.get<std::string>("PFParticleModule", "pandora");
    m_trackfitLabel = pset.get<std::string>("TrackFitModule", "trackfit");
    m_hitfinderLabel = pset.get<std::string>("HitFinderModule", "gaushit");
    m_geantModuleLabel = pset.get<std::string>("GeantModule", "largeant");

    m_useDaughterPFParticles = pset.get<bool>("UseDaughterPFParticles", true);
    m_useDaughterMCParticles = pset.get<bool>("UseDaughterMCParticles", true);

    m_cosmicContainmentCut = pset.get<double>("CosmicContainmentCut", 5.0);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void PFParticleCosmicAna::beginJob()
  {
    mf::LogDebug("LArPandora") << " *** PFParticleCosmicAna::beginJob() *** " << std::endl;

    //
    art::ServiceHandle<art::TFileService const> tfs;

    m_pRecoTree = tfs->make<TTree>("recoTree", "LAr Cosmic Reco Tree");
    m_pRecoTree->Branch("run", &m_run, "run/I");
    m_pRecoTree->Branch("event", &m_event, "event/I");
    m_pRecoTree->Branch("index", &m_index, "index/I");
    m_pRecoTree->Branch("self", &m_self, "self/I");
    m_pRecoTree->Branch("pdgCode", &m_pdgCode, "pdgCode/I");
    m_pRecoTree->Branch("isTrackLike", &m_isTrackLike, "isTrackLike/I");
    m_pRecoTree->Branch("isPrimary", &m_isPrimary, "isPrimary/I");
    m_pRecoTree->Branch("cosmicScore", &m_cosmicScore, "cosmicScore/F");
    m_pRecoTree->Branch("trackVtxX", &m_trackVtxX, "trackVtxX/F");
    m_pRecoTree->Branch("trackVtxY", &m_trackVtxY, "trackVtxY/F");
    m_pRecoTree->Branch("trackVtxZ", &m_trackVtxZ, "trackVtxZ/F");
    m_pRecoTree->Branch("trackEndX", &m_trackEndX, "trackEndX/F");
    m_pRecoTree->Branch("trackEndY", &m_trackEndY, "trackEndY/F");
    m_pRecoTree->Branch("trackEndZ", &m_trackEndZ, "trackEndZ/F");
    m_pRecoTree->Branch("trackVtxDirX", &m_trackVtxDirX, "trackVtxDirX/F");
    m_pRecoTree->Branch("trackVtxDirY", &m_trackVtxDirY, "trackVtxDirY/F");
    m_pRecoTree->Branch("trackVtxDirZ", &m_trackVtxDirZ, "trackVtxDirZ/F");
    m_pRecoTree->Branch("trackEndDirX", &m_trackEndDirX, "trackEndDirX/F");
    m_pRecoTree->Branch("trackEndDirY", &m_trackEndDirY, "trackEndDirY/F");
    m_pRecoTree->Branch("trackEndDirZ", &m_trackEndDirZ, "trackEndDirZ/F");
    m_pRecoTree->Branch("trackLength", &m_trackLength, "trackLength/F");
    m_pRecoTree->Branch("trackWidthX", &m_trackWidthX, "trackWidthX/F");
    m_pRecoTree->Branch("trackWidthY", &m_trackWidthY, "trackWidthY/F");
    m_pRecoTree->Branch("trackWidthZ", &m_trackWidthZ, "trackWidthZ/F");
    m_pRecoTree->Branch("trackVtxDeltaYZ", &m_trackVtxDeltaYZ, "trackVtxDeltaYZ/F");
    m_pRecoTree->Branch("trackEndDeltaYZ", &m_trackEndDeltaYZ, "trackEndDeltaYZ/F");
    m_pRecoTree->Branch("trackVtxContained", &m_trackVtxContained, "trackVtxContained/I");
    m_pRecoTree->Branch("trackEndContained", &m_trackEndContained, "trackEndContained/I");
    m_pRecoTree->Branch("nTracks", &m_nTracks, "nTracks/I");
    m_pRecoTree->Branch("nHits", &m_nHits, "nHits/I");

    m_pTrueTree = tfs->make<TTree>("trueTree", "LAr Cosmic True Tree");
    m_pTrueTree->Branch("run", &m_run, "run/I");
    m_pTrueTree->Branch("event", &m_event, "event/I");
    m_pTrueTree->Branch("nHits", &m_nHits, "nHits/I");
    m_pTrueTree->Branch("nNeutrinoHits", &m_nNeutrinoHits, "nNeutrinoHits/I");
    m_pTrueTree->Branch(
      "nNeutrinoHitsFullyTagged", &m_nNeutrinoHitsFullyTagged, "nNeutrinoHitsFullyTagged/I");
    m_pTrueTree->Branch(
      "nNeutrinoHitsSemiTagged", &m_nNeutrinoHitsSemiTagged, "nNeutrinoHitsSemiTagged/I");
    m_pTrueTree->Branch(
      "nNeutrinoHitsNotTagged", &m_nNeutrinoHitsNotTagged, "nNeutrinoHitsNotTagged/I");
    m_pTrueTree->Branch("nNeutrinoHitsNotReconstructed",
                        &m_nNeutrinoHitsNotReconstructed,
                        "nNeutrinoHitsNotReconstructed/I");
    m_pTrueTree->Branch(
      "nNeutrinoHitsReconstructed", &m_nNeutrinoHitsReconstructed, "nNeutrinoHitsReconstructed/I");
    m_pTrueTree->Branch("nCosmicHits", &m_nCosmicHits, "nCosmicHits/I");
    m_pTrueTree->Branch(
      "nCosmicHitsFullyTagged", &m_nCosmicHitsFullyTagged, "nCosmicHitsFullyTagged/I");
    m_pTrueTree->Branch(
      "nCosmicHitsSemiTagged", &m_nCosmicHitsSemiTagged, "nCosmicHitsSemiTagged/I");
    m_pTrueTree->Branch("nCosmicHitsNotTagged", &m_nCosmicHitsNotTagged, "nCosmicHitsNotTagged/I");
    m_pTrueTree->Branch("nCosmicHitsNotReconstructed",
                        &m_nCosmicHitsNotReconstructed,
                        "nCosmicHitsNotReconstructed/I");
    m_pTrueTree->Branch(
      "nCosmicHitsReconstructed", &m_nCosmicHitsReconstructed, "nCosmicHitsReconstructed/I");
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void PFParticleCosmicAna::endJob() {}

  //------------------------------------------------------------------------------------------------------------------------------------------

  void PFParticleCosmicAna::analyze(const art::Event& evt)
  {
    std::cout << " *** PFParticleCosmicAna::analyze(...) *** " << std::endl;

    //
    // Note: I've made this is MicroBooNE-only module
    //

    m_run = evt.run();
    m_event = evt.id().event();

    std::cout << "  Run: " << m_run << std::endl;
    std::cout << "  Event: " << m_event << std::endl;

    // Collect True Particles
    // ======================
    HitVector hitVector;
    MCTruthToMCParticles truthToParticles;
    MCParticlesToMCTruth particlesToTruth;
    MCParticlesToHits trueParticlesToHits;
    HitsToMCParticles trueHitsToParticles;

    LArPandoraHelper::CollectHits(evt, m_hitfinderLabel, hitVector);
    LArPandoraHelper::CollectMCParticles(
      evt, m_geantModuleLabel, truthToParticles, particlesToTruth);
    LArPandoraHelper::BuildMCParticleHitMaps(evt,
                                             m_geantModuleLabel,
                                             hitVector,
                                             trueParticlesToHits,
                                             trueHitsToParticles,
                                             (m_useDaughterMCParticles ?
                                                LArPandoraHelper::kAddDaughters :
                                                LArPandoraHelper::kIgnoreDaughters));

    // Collect Reco Particles
    // ======================
    PFParticleVector recoParticleVector;
    PFParticlesToHits recoParticlesToHits;
    HitsToPFParticles recoHitsToParticles;

    LArPandoraHelper::CollectPFParticles(evt, m_particleLabel, recoParticleVector);
    LArPandoraHelper::BuildPFParticleHitMaps(evt,
                                             m_particleLabel,
                                             m_particleLabel,
                                             recoParticlesToHits,
                                             recoHitsToParticles,
                                             (m_useDaughterPFParticles ?
                                                LArPandoraHelper::kAddDaughters :
                                                LArPandoraHelper::kIgnoreDaughters));

    std::cout << "  PFParticles: " << recoParticleVector.size() << std::endl;

    // Collect Reco Tracks
    // ===================
    TrackVector recoTrackVector;
    PFParticlesToTracks recoParticlesToTracks;
    LArPandoraHelper::CollectTracks(evt, m_trackfitLabel, recoTrackVector, recoParticlesToTracks);

    // Collect Cosmic Tags
    // =====================
    CosmicTagVector recoCosmicTagVector;
    TracksToCosmicTags recoTracksToCosmicTags;
    LArPandoraHelper::CollectCosmicTags(
      evt, m_cosmicLabel, recoCosmicTagVector, recoTracksToCosmicTags);

    // Analyse Reconstructed Particles
    // ===============================
    this->FillRecoTree(recoParticlesToHits, recoParticlesToTracks, recoTracksToCosmicTags);

    // Analyse True Hits
    // =================
    this->FillTrueTree(hitVector,
                       trueHitsToParticles,
                       recoHitsToParticles,
                       particlesToTruth,
                       recoParticlesToTracks,
                       recoTracksToCosmicTags);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void PFParticleCosmicAna::FillRecoTree(const PFParticlesToHits& recoParticlesToHits,
                                         const PFParticlesToTracks& recoParticlesToTracks,
                                         const TracksToCosmicTags& recoTracksToCosmicTags)
  {
    // Set up Geometry Service
    // =======================
    auto const& tpc = art::ServiceHandle<geo::Geometry const>()->TPC({0, 0});

    const double xmin(0.0);
    const double xmax(2.0 * tpc.HalfWidth());
    const double ymin(-tpc.HalfHeight());
    const double ymax(+tpc.HalfHeight());
    const double zmin(0.0);
    const double zmax(tpc.Length());
    const double xyzCut(m_cosmicContainmentCut);

    m_index = 0;

    m_self = 0;
    m_pdgCode = 0;
    m_isTrackLike = 0;
    m_isPrimary = 0;
    m_cosmicScore = 0.f;

    m_trackVtxX = 0.f;
    m_trackVtxY = 0.f;
    m_trackVtxZ = 0.f;
    m_trackEndX = 0.f;
    m_trackEndY = 0.f;
    m_trackEndZ = 0.f;
    m_trackVtxDirX = 0.f;
    m_trackVtxDirY = 0.f;
    m_trackVtxDirZ = 0.f;
    m_trackEndDirX = 0.f;
    m_trackEndDirY = 0.f;
    m_trackEndDirZ = 0.f;
    m_trackLength = 0.f;
    m_trackWidthX = 0.f;
    m_trackWidthY = 0.f;
    m_trackWidthZ = 0.f;
    m_trackVtxDeltaYZ = 0.f;
    m_trackEndDeltaYZ = 0.f;

    m_trackVtxContained = 0;
    m_trackEndContained = 0;

    m_nTracks = 0;
    m_nHits = 0;

    // Loop over Reco Particles
    // ========================
    for (PFParticlesToHits::const_iterator iter1 = recoParticlesToHits.begin(),
                                           iterEnd1 = recoParticlesToHits.end();
         iter1 != iterEnd1;
         ++iter1) {
      const art::Ptr<recob::PFParticle> recoParticle = iter1->first;

      const HitVector& hitVector = iter1->second;
      if (hitVector.empty()) continue;

      PFParticlesToTracks::const_iterator iter2 = recoParticlesToTracks.find(recoParticle);
      if (recoParticlesToTracks.end() == iter2) continue;

      const TrackVector& trackVector = iter2->second;
      if (trackVector.empty()) continue;

      m_nHits = hitVector.size();
      m_nTracks = trackVector.size();

      m_self = recoParticle->Self();
      m_pdgCode = recoParticle->PdgCode();
      m_isPrimary = recoParticle->IsPrimary();
      m_isTrackLike = LArPandoraHelper::IsTrack(recoParticle);
      m_cosmicScore =
        this->GetCosmicScore(recoParticle, recoParticlesToTracks, recoTracksToCosmicTags);

      m_trackVtxX = 0.f;
      m_trackVtxY = 0.f;
      m_trackVtxZ = 0.f;
      m_trackEndX = 0.f;
      m_trackEndY = 0.f;
      m_trackEndZ = 0.f;
      m_trackVtxDirX = 0.f;
      m_trackVtxDirY = 0.f;
      m_trackVtxDirZ = 0.f;
      m_trackEndDirX = 0.f;
      m_trackEndDirY = 0.f;
      m_trackEndDirZ = 0.f;
      m_trackLength = 0.f;
      m_trackWidthX = 0.f;
      m_trackWidthY = 0.f;
      m_trackWidthZ = 0.f;
      m_trackVtxDeltaYZ = 0.f;
      m_trackEndDeltaYZ = 0.f;

      m_trackVtxContained = 0;
      m_trackEndContained = 0;

      for (TrackVector::const_iterator iter3 = trackVector.begin(), iterEnd3 = trackVector.end();
           iter3 != iterEnd3;
           ++iter3) {
        const art::Ptr<recob::Track> track = *iter3;
        const float trackLength(track->Length());

        if (trackLength < m_trackLength) continue;

        m_trackLength = trackLength;

        const auto& trackVtxPosition = track->Vertex();
        const auto& trackVtxDirection = track->VertexDirection();
        const auto& trackEndPosition = track->End();
        const auto& trackEndDirection = track->EndDirection();

        m_trackVtxX = trackVtxPosition.x();
        m_trackVtxY = trackVtxPosition.y();
        m_trackVtxZ = trackVtxPosition.z();
        m_trackVtxDirX = trackVtxDirection.x();
        m_trackVtxDirY = trackVtxDirection.y();
        m_trackVtxDirZ = trackVtxDirection.z();
        m_trackEndX = trackEndPosition.x();
        m_trackEndY = trackEndPosition.y();
        m_trackEndZ = trackEndPosition.z();
        m_trackEndDirX = trackEndDirection.x();
        m_trackEndDirY = trackEndDirection.y();
        m_trackEndDirZ = trackEndDirection.z();

        m_trackWidthX = std::fabs(m_trackEndX - m_trackVtxX);
        m_trackWidthY = std::fabs(m_trackEndY - m_trackVtxY);
        m_trackWidthZ = std::fabs(m_trackEndZ - m_trackVtxZ);

        m_trackVtxDeltaYZ =
          std::min((ymax - m_trackVtxY), std::min((m_trackVtxZ - zmin), (zmax - m_trackVtxZ)));
        m_trackEndDeltaYZ =
          std::min((m_trackEndY - ymin), std::min((m_trackEndZ - zmin), (zmax - m_trackEndZ)));

        m_trackVtxContained = ((m_trackVtxX > xmin + xyzCut && m_trackVtxX < xmax - xyzCut) &&
                               (m_trackVtxY > ymin + xyzCut && m_trackVtxY < ymax - xyzCut) &&
                               (m_trackVtxZ > zmin + xyzCut && m_trackVtxZ < zmax - xyzCut));
        m_trackEndContained = ((m_trackEndX > xmin + xyzCut && m_trackEndX < xmax - xyzCut) &&
                               (m_trackEndY > ymin + xyzCut && m_trackEndY < ymax - xyzCut) &&
                               (m_trackEndZ > zmin + xyzCut && m_trackEndZ < zmax - xyzCut));
      }

      std::cout << "   PFParticle: [" << m_index << "] nHits=" << m_nHits
                << ", nTracks=" << m_nTracks << ", cosmicScore=" << m_cosmicScore << std::endl;

      m_pRecoTree->Fill();
      ++m_index;
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void PFParticleCosmicAna::FillTrueTree(const HitVector& hitVector,
                                         const HitsToMCParticles& trueHitsToParticles,
                                         const HitsToPFParticles& recoHitsToParticles,
                                         const MCParticlesToMCTruth& particlesToTruth,
                                         const PFParticlesToTracks& particlesToTracks,
                                         const TracksToCosmicTags& tracksToCosmicTags)
  {
    m_nHits = 0;

    m_nNeutrinoHits = 0;
    m_nNeutrinoHitsFullyTagged = 0;
    m_nNeutrinoHitsSemiTagged = 0;
    m_nNeutrinoHitsNotTagged = 0;
    m_nNeutrinoHitsNotReconstructed = 0;
    m_nNeutrinoHitsReconstructed = 0;

    m_nCosmicHits = 0;
    m_nCosmicHitsFullyTagged = 0;
    m_nCosmicHitsSemiTagged = 0;
    m_nCosmicHitsNotTagged = 0;
    m_nCosmicHitsNotReconstructed = 0;
    m_nCosmicHitsReconstructed = 0;

    for (HitVector::const_iterator iter2 = hitVector.begin(), iterEnd2 = hitVector.end();
         iter2 != iterEnd2;
         ++iter2) {
      const art::Ptr<recob::Hit> hit = *iter2;

      HitsToMCParticles::const_iterator iter3 = trueHitsToParticles.find(hit);
      if (trueHitsToParticles.end() == iter3) continue;

      const art::Ptr<simb::MCParticle> trueParticle = iter3->second;

      MCParticlesToMCTruth::const_iterator iter4 = particlesToTruth.find(trueParticle);
      if (particlesToTruth.end() == iter4)
        throw cet::exception("LArPandora") << " PFParticleCosmicAna::analyze --- Found a true "
                                              "particle without any ancestry information ";

      const art::Ptr<simb::MCTruth> truth = iter4->second;

      float cosmicScore(-0.2);

      HitsToPFParticles::const_iterator iter5 = recoHitsToParticles.find(hit);
      if (recoHitsToParticles.end() != iter5) {
        const art::Ptr<recob::PFParticle> particle = iter5->second;
        cosmicScore = this->GetCosmicScore(particle, particlesToTracks, tracksToCosmicTags);
      }

      ++m_nHits;

      if (truth->NeutrinoSet()) {
        ++m_nNeutrinoHits;

        if (cosmicScore >= 0)
          ++m_nNeutrinoHitsReconstructed;
        else
          ++m_nNeutrinoHitsNotReconstructed;

        if (cosmicScore > 0.51)
          ++m_nNeutrinoHitsFullyTagged;
        else if (cosmicScore > 0.39)
          ++m_nNeutrinoHitsSemiTagged;
        else
          ++m_nNeutrinoHitsNotTagged;
      }
      else {
        ++m_nCosmicHits;

        if (cosmicScore >= 0)
          ++m_nCosmicHitsReconstructed;
        else
          ++m_nCosmicHitsNotReconstructed;

        if (cosmicScore > 0.51)
          ++m_nCosmicHitsFullyTagged;
        else if (cosmicScore > 0.39)
          ++m_nCosmicHitsSemiTagged;
        else
          ++m_nCosmicHitsNotTagged;
      }
    }

    m_pTrueTree->Fill();
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  float PFParticleCosmicAna::GetCosmicScore(const art::Ptr<recob::PFParticle> particle,
                                            const PFParticlesToTracks& recoParticlesToTracks,
                                            const TracksToCosmicTags& recoTracksToCosmicTags) const
  {
    float cosmicScore(0.f);

    // Get cosmic tags associated with this particle
    PFParticlesToTracks::const_iterator iter2 = recoParticlesToTracks.find(particle);
    if (recoParticlesToTracks.end() != iter2) {
      for (TrackVector::const_iterator iter3 = iter2->second.begin(),
                                       iterEnd3 = iter2->second.end();
           iter3 != iterEnd3;
           ++iter3) {
        const art::Ptr<recob::Track> track = *iter3;

        TracksToCosmicTags::const_iterator iter4 = recoTracksToCosmicTags.find(track);
        if (recoTracksToCosmicTags.end() != iter4) {
          for (CosmicTagVector::const_iterator iter5 = iter4->second.begin(),
                                               iterEnd5 = iter4->second.end();
               iter5 != iterEnd5;
               ++iter5) {
            const art::Ptr<anab::CosmicTag> cosmicTag = *iter5;
            if (cosmicTag->CosmicScore() > cosmicScore) cosmicScore = cosmicTag->CosmicScore();
          }
        }
      }
    }

    return cosmicScore;
  }

} //namespace lar_pandora
