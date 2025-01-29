/**
 *  @file   larpandora/LArPandoraAnalysis/PFParticleHitDumper_module.cc
 *
 *  @brief  Analysis module for created particles
 *
 */

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "TTree.h"

#include "larcore/Geometry/WireReadout.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora {

  /**
   *  @brief  PFParticleHitDumper class
   */
  class PFParticleHitDumper : public art::EDAnalyzer {
  public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    PFParticleHitDumper(fhicl::ParameterSet const& pset);

    /**
     *  @brief  Destructor
     */
    virtual ~PFParticleHitDumper();

    void beginJob();
    void endJob();
    void analyze(const art::Event& evt);
    void reconfigure(fhicl::ParameterSet const& pset);

  private:
    /**
     *  @brief Store 3D track hits
     *
     *  @param particlesToTracks mapping between 3D track hits and PFParticles
     */
    void FillRecoTracks(const PFParticlesToTracks& particlesToTracks);

    /**
     *  @brief Store 3D hits
     *
     *  @param particleVector the input vector of PFParticles
     *  @param particlesToSpacePoints mapping between 3D hits and PFParticles
     *  @param spacePointsToHits mapping between 3D hits and 2D hits
     */
    void FillReco3D(const PFParticleVector& particleVector,
                    const PFParticlesToSpacePoints& particlesToSpacePoints,
                    const SpacePointsToHits& spacePointsToHits);

    /**
     *  @brief Store 2D hits
     *
     *  @param hitVector the input vector of 2D hits
     *  @param hitsToParticles mapping between 2D hits and PFParticles
     */
    void FillReco2D(const art::Event& event,
                    const HitVector& hitVector,
                    const HitsToPFParticles& hitsToParticles);

    /**
     *  @brief Store number of 2D hits associated to PFParticle in different ways
     *
     *  @param particleVector the input vector of PFParticles
     *  @param particlesToHits mapping between PFParticles and 2D hits through space points
     *  @param particlesToHitsClusters mapping between PFParticles and 2D hits through clusters
     *  @param particlesToTracks mapping between PFParticles and tracks
     *  @param particlesToShowers mapping between PFParticles and showers
     */
    void FillAssociated2DHits(const art::Event& evt,
                              const PFParticleVector& particleVector,
                              const PFParticlesToHits& particlesToHits,
                              const PFParticlesToHits& particlesToHitsClusters,
                              const PFParticlesToTracks& particlesToTracks,
                              const TracksToHits& tracksToHits,
                              const PFParticlesToShowers& particlesToShowers,
                              const ShowersToHits& showersToHits);

    /**
     *  @brief Store raw data
     *
     *  @param wireVector the input vector of reconstructed wires
     */
    void FillRecoWires(const art::Event& event, const WireVector& wireVector);

    /**
     *  @brief Conversion from wire ID to U/V/W coordinate
     *
     *  @param wireID the input wire ID
     */
    double GetUVW(const geo::WireID& wireID) const;

    /**
     *  @brief Convert from (Y,Z) to U coordinate
     *
     *  @param cstat the cryostat
     *  @param tpc the tpc
     *  @param y the y coordinate
     *  @param z the z coordinate
     */
    double YZtoU(const unsigned int cstat,
                 const unsigned int tpc,
                 const double y,
                 const double z) const;

    /**
     *  @brief Convert from (Y,Z) to V coordinate
     *
     *  @param cstat the crystat
     *  @param tpc the tpc
     *  @param y the y coordinate
     *  @param z the z coordinate
     */
    double YZtoV(const unsigned int cstat,
                 const unsigned int tpc,
                 const double y,
                 const double z) const;

    TTree* m_pRecoTracks;     ///<
    TTree* m_pReco3D;         ///<
    TTree* m_pReco2D;         ///<
    TTree* m_pRecoComparison; ///<
    TTree* m_pRecoWire;       ///<

    int m_run;      ///<
    int m_event;    ///<
    int m_particle; ///<
    int m_primary;  ///<
    int m_pdgcode;  ///<

    int m_cstat; ///<
    int m_tpc;   ///<
    int m_plane; ///<
    int m_wire;  ///<

    double m_u; ///<
    double m_v; ///<
    double m_w; ///<
    double m_x; ///<
    double m_y; ///<
    double m_z; ///<
    double m_q; ///<

    int m_hitsFromSpacePoints;   ///<
    int m_hitsFromClusters;      ///<
    int m_hitsFromTrackOrShower; ///<

    std::string m_calwireLabel;    ///<
    std::string m_hitfinderLabel;  ///<
    std::string m_clusterLabel;    ///<
    std::string m_spacepointLabel; ///<
    std::string m_particleLabel;   ///<
    std::string m_trackLabel;      ///<
    std::string m_showerLabel;     ///<

    bool m_storeWires; ///<
    bool m_printDebug; ///< switch for print statements (TODO: use message service!)
  };

  DEFINE_ART_MODULE(PFParticleHitDumper)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Wire.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include <iostream>

namespace lar_pandora {

  PFParticleHitDumper::PFParticleHitDumper(fhicl::ParameterSet const& pset) : art::EDAnalyzer(pset)
  {
    this->reconfigure(pset);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  PFParticleHitDumper::~PFParticleHitDumper() {}

  //------------------------------------------------------------------------------------------------------------------------------------------

  void PFParticleHitDumper::reconfigure(fhicl::ParameterSet const& pset)
  {
    m_storeWires = pset.get<bool>("StoreWires", false);
    m_trackLabel = pset.get<std::string>("TrackModule", "pandoraTrack");
    m_showerLabel = pset.get<std::string>("ShowerModule", "pandoraShower");
    m_particleLabel = pset.get<std::string>("PFParticleModule", "pandora");
    m_spacepointLabel = pset.get<std::string>("SpacePointModule", "pandora");
    m_clusterLabel = pset.get<std::string>("ClusterModule", "pandora");
    m_hitfinderLabel = pset.get<std::string>("HitFinderModule", "gaushit");
    m_calwireLabel = pset.get<std::string>("CalWireModule", "caldata");
    m_printDebug = pset.get<bool>("PrintDebug", false);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void PFParticleHitDumper::beginJob()
  {
    mf::LogDebug("LArPandora") << " *** PFParticleHitDumper::beginJob() *** " << std::endl;

    //
    art::ServiceHandle<art::TFileService const> tfs;

    m_pRecoTracks = tfs->make<TTree>("pandoraTracks", "LAr Reco Tracks");
    m_pRecoTracks->Branch("run", &m_run, "run/I");
    m_pRecoTracks->Branch("event", &m_event, "event/I");
    m_pRecoTracks->Branch("particle", &m_particle, "particle/I");
    m_pRecoTracks->Branch("x", &m_x, "x/D");
    m_pRecoTracks->Branch("y", &m_y, "y/D");
    m_pRecoTracks->Branch("z", &m_z, "z/D");

    m_pReco3D = tfs->make<TTree>("pandora3D", "LAr Reco 3D");
    m_pReco3D->Branch("run", &m_run, "run/I");
    m_pReco3D->Branch("event", &m_event, "event/I");
    m_pReco3D->Branch("particle", &m_particle, "particle/I");
    m_pReco3D->Branch("primary", &m_primary, "primary/I");
    m_pReco3D->Branch("pdgcode", &m_pdgcode, "pdgcode/I");
    m_pReco3D->Branch("cstat", &m_cstat, "cstat/I");
    m_pReco3D->Branch("tpc", &m_tpc, "tpc/I");
    m_pReco3D->Branch("plane", &m_plane, "plane/I");
    m_pReco3D->Branch("x", &m_x, "x/D");
    m_pReco3D->Branch("y", &m_y, "y/D");
    m_pReco3D->Branch("u", &m_u, "u/D");
    m_pReco3D->Branch("v", &m_v, "v/D");
    m_pReco3D->Branch("z", &m_z, "z/D");

    m_pReco2D = tfs->make<TTree>("pandora2D", "LAr Reco 2D");
    m_pReco2D->Branch("run", &m_run, "run/I");
    m_pReco2D->Branch("event", &m_event, "event/I");
    m_pReco2D->Branch("particle", &m_particle, "particle/I");
    m_pReco2D->Branch("pdgcode", &m_pdgcode, "pdgcode/I");
    m_pReco2D->Branch("cstat", &m_cstat, "cstat/I");
    m_pReco2D->Branch("tpc", &m_tpc, "tpc/I");
    m_pReco2D->Branch("plane", &m_plane, "plane/I");
    m_pReco2D->Branch("wire", &m_wire, "wire/I");
    m_pReco2D->Branch("x", &m_x, "x/D");
    m_pReco2D->Branch("w", &m_w, "w/D");
    m_pReco2D->Branch("q", &m_q, "q/D");

    m_pRecoComparison = tfs->make<TTree>("pandora2Dcomparison", "LAr Reco 2D (comparison)");
    m_pRecoComparison->Branch("run", &m_run, "run/I");
    m_pRecoComparison->Branch("event", &m_event, "event/I");
    m_pRecoComparison->Branch("particle", &m_particle, "particle/I");
    m_pRecoComparison->Branch("pdgcode", &m_pdgcode, "pdgcode/I");
    m_pRecoComparison->Branch(
      "hitsFromSpacePoints", &m_hitsFromSpacePoints, "hitsFromSpacePoints/I");
    m_pRecoComparison->Branch("hitsFromClusters", &m_hitsFromClusters, "hitsFromClusters/I");
    m_pRecoComparison->Branch(
      "hitsFromTrackOrShower", &m_hitsFromTrackOrShower, "hitsFromTrackOrShower/I");

    m_pRecoWire = tfs->make<TTree>("rawdata", "LAr Reco Wires");
    m_pRecoWire->Branch("run", &m_run, "run/I");
    m_pRecoWire->Branch("event", &m_event, "event/I");
    m_pRecoWire->Branch("cstat", &m_cstat, "cstat/I");
    m_pRecoWire->Branch("tpc", &m_tpc, "tpc/I");
    m_pRecoWire->Branch("plane", &m_plane, "plane/I");
    m_pRecoWire->Branch("wire", &m_wire, "wire/I");
    m_pRecoWire->Branch("x", &m_x, "x/D");
    m_pRecoWire->Branch("w", &m_w, "w/D");
    m_pRecoWire->Branch("q", &m_q, "q/D");
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void PFParticleHitDumper::endJob() {}

  //------------------------------------------------------------------------------------------------------------------------------------------

  void PFParticleHitDumper::analyze(const art::Event& evt)
  {
    if (m_printDebug) std::cout << " *** PFParticleHitDumper::analyze(...) *** " << std::endl;

    m_run = evt.run();
    m_event = evt.id().event();

    m_particle = -1;
    m_primary = 0;
    m_pdgcode = 0;

    m_cstat = 0;
    m_tpc = 0;
    m_plane = 0;
    m_wire = 0;

    m_x = 0.0;
    m_y = 0.0;
    m_u = 0.0;
    m_v = 0.0;
    m_z = 0.0;
    m_w = 0.0;
    m_q = 0.0;

    if (m_printDebug) {
      std::cout << "  Run: " << m_run << std::endl;
      std::cout << "  Event: " << m_event << std::endl;
    }

    // Get particles, tracks, space points, hits (and wires)
    // ====================================================
    TrackVector trackVector, trackVectorExtra;
    ShowerVector showerVector, showerVectorExtra;
    PFParticleVector particleVector;
    SpacePointVector spacePointVector;
    HitVector hitVector;
    WireVector wireVector;

    PFParticlesToTracks particlesToTracks;
    PFParticlesToShowers particlesToShowers;
    PFParticlesToSpacePoints particlesToSpacePoints;
    PFParticlesToHits particlesToHits, particlesToHitsClusters;
    TracksToHits tracksToHits;
    ShowersToHits showersToHits;
    HitsToPFParticles hitsToParticles, hitsToParticlesClusters;
    SpacePointsToHits spacePointsToHits;

    LArPandoraHelper::CollectHits(evt, m_hitfinderLabel, hitVector);
    LArPandoraHelper::CollectSpacePoints(
      evt, m_spacepointLabel, spacePointVector, spacePointsToHits);
    LArPandoraHelper::CollectTracks(evt, m_trackLabel, trackVector, particlesToTracks);
    LArPandoraHelper::CollectTracks(evt, m_trackLabel, trackVectorExtra, tracksToHits);
    LArPandoraHelper::CollectShowers(evt, m_showerLabel, showerVector, particlesToShowers);
    LArPandoraHelper::CollectShowers(evt, m_showerLabel, showerVectorExtra, showersToHits);
    LArPandoraHelper::CollectPFParticles(
      evt, m_particleLabel, particleVector, particlesToSpacePoints);
    LArPandoraHelper::BuildPFParticleHitMaps(evt,
                                             m_particleLabel,
                                             m_spacepointLabel,
                                             particlesToHits,
                                             hitsToParticles,
                                             LArPandoraHelper::DaughterMode::kUseDaughters,
                                             false);
    LArPandoraHelper::BuildPFParticleHitMaps(
      evt, m_particleLabel, m_clusterLabel, particlesToHitsClusters, hitsToParticlesClusters);

    if (m_storeWires) LArPandoraHelper::CollectWires(evt, m_calwireLabel, wireVector);

    if (m_printDebug) std::cout << "  PFParticles: " << particleVector.size() << std::endl;

    // Loop over Tracks (Fill 3D Track Tree)
    // =====================================
    if (m_printDebug) std::cout << "   PFParticleHitDumper::FillRecoTracks(...) " << std::endl;
    this->FillRecoTracks(particlesToTracks);

    // Loop over PFParticles (Fill 3D Reco Tree)
    // =========================================
    if (m_printDebug) std::cout << "   PFParticleHitDumper::FillReco3D(...) " << std::endl;
    this->FillReco3D(particleVector, particlesToSpacePoints, spacePointsToHits);

    // Loop over Hits (Fill 2D Reco Tree)
    // ==================================
    if (m_printDebug) std::cout << "   PFParticleHitDumper::FillReco2D(...) " << std::endl;
    this->FillReco2D(evt, hitVector, hitsToParticles);

    // Loop over Hits (Fill Associated 2D Hits Tree)
    // =============================================
    if (m_printDebug)
      std::cout << "   PFParticleHitDumper::FillAssociated2DHits(...) " << std::endl;
    this->FillAssociated2DHits(evt,
                               particleVector,
                               particlesToHits,
                               particlesToHitsClusters,
                               particlesToTracks,
                               tracksToHits,
                               particlesToShowers,
                               showersToHits);

    // Loop over Wires (Fill Reco Wire Tree)
    // =====================================
    if (m_printDebug) std::cout << "   PFParticleHitDumper::FillRecoWires(...) " << std::endl;
    this->FillRecoWires(evt, wireVector);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void PFParticleHitDumper::FillRecoTracks(const PFParticlesToTracks& particlesToTracks)
  {
    // Initialise variables
    m_particle = -1;
    m_x = 0.0;
    m_y = 0.0;
    m_z = 0.0;

    // Create dummy entry if there are no particles
    if (particlesToTracks.empty()) { m_pRecoTracks->Fill(); }

    // Loop over tracks
    for (PFParticlesToTracks::const_iterator iter = particlesToTracks.begin(),
                                             iterEnd = particlesToTracks.end();
         iter != iterEnd;
         ++iter) {
      const art::Ptr<recob::PFParticle> particle = iter->first;
      const TrackVector& trackVector = iter->second;

      m_particle = particle->Self();

      if (!trackVector.empty()) {
        if (trackVector.size() != 1 && m_printDebug)
          std::cout << " Warning: Found particle with more than one associated track " << std::endl;

        const art::Ptr<recob::Track> track = *(trackVector.begin());

        if (m_printDebug)
          std::cout << "    PFPARTICLE [" << m_particle << "] (" << track->NumberTrajectoryPoints()
                    << " Trajectory Points)" << std::endl;

        for (unsigned int p = 0; p < track->NumberTrajectoryPoints(); ++p) {
          const auto position(track->LocationAtPoint(p));
          m_x = position.x();
          m_y = position.y();
          m_z = position.z();

          m_pRecoTracks->Fill();
        }
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void PFParticleHitDumper::FillReco3D(const PFParticleVector& particleVector,
                                       const PFParticlesToSpacePoints& particlesToSpacePoints,
                                       const SpacePointsToHits& spacePointsToHits)
  {
    // Initialise variables
    m_particle = -1;
    m_primary = 0;
    m_pdgcode = 0;
    m_cstat = 0;
    m_tpc = 0;
    m_plane = 0;
    m_x = 0.0;
    m_u = 0.0;
    m_v = 0.0;
    m_y = 0.0;
    m_z = 0.0;

    // Create dummy entry if there are no particles
    if (particleVector.empty()) { m_pReco3D->Fill(); }

    // Store associations between particle and particle ID
    PFParticleMap theParticleMap;

    for (unsigned int i = 0; i < particleVector.size(); ++i) {
      const art::Ptr<recob::PFParticle> particle = particleVector.at(i);
      theParticleMap[particle->Self()] = particle;
    }

    // Loop over particles
    for (PFParticlesToSpacePoints::const_iterator iter1 = particlesToSpacePoints.begin(),
                                                  iterEnd1 = particlesToSpacePoints.end();
         iter1 != iterEnd1;
         ++iter1) {
      const art::Ptr<recob::PFParticle> particle = iter1->first;
      const SpacePointVector& spacepoints = iter1->second;

      m_particle = particle->Self();
      m_pdgcode = particle->PdgCode();
      m_primary = 0;

      if (particle->IsPrimary()) { m_primary = 1; }
      else {
        const size_t parentID(particle->Parent());
        PFParticleMap::const_iterator pIter = theParticleMap.find(parentID);
        if (theParticleMap.end() == pIter)
          throw cet::exception("LArPandora")
            << " PFParticleHitDumper::analyze --- Found particle with ID code";

        const art::Ptr<recob::PFParticle> particleParent = pIter->second;
        if (LArPandoraHelper::IsNeutrino(particleParent)) m_primary = 1;
      }

      if (m_printDebug)
        std::cout << "    PFPARTICLE [" << m_particle << "] [Primary=" << m_primary << "] ("
                  << spacepoints.size() << " Space Points)" << std::endl;

      for (unsigned int j = 0; j < spacepoints.size(); ++j) {
        const art::Ptr<recob::SpacePoint> spacepoint = spacepoints.at(j);

        m_x = spacepoint->XYZ()[0];
        m_y = spacepoint->XYZ()[1];
        m_z = spacepoint->XYZ()[2];

        SpacePointsToHits::const_iterator iter2 = spacePointsToHits.find(spacepoint);
        if (spacePointsToHits.end() == iter2)
          throw cet::exception("LArPandora")
            << " PFParticleHitDumper::analyze --- Found space point without associated hit";

        const art::Ptr<recob::Hit> hit = iter2->second;
        const geo::WireID& wireID(hit->WireID());

        m_cstat = wireID.Cryostat;
        m_tpc = wireID.TPC;
        m_plane = wireID.Plane;

        m_u = this->YZtoU(m_cstat, m_tpc, m_y, m_z);
        m_v = this->YZtoV(m_cstat, m_tpc, m_y, m_z);

        m_pReco3D->Fill();
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void PFParticleHitDumper::FillAssociated2DHits(const art::Event& evt,
                                                 const PFParticleVector& particleVector,
                                                 const PFParticlesToHits& particlesToHits,
                                                 const PFParticlesToHits& particlesToHitsClusters,
                                                 const PFParticlesToTracks& particlesToTracks,
                                                 const TracksToHits& tracksToHits,
                                                 const PFParticlesToShowers& particlesToShowers,
                                                 const ShowersToHits& showersToHits)
  {
    // Create dummy entry if there are no 2D hits
    if (particleVector.empty()) { m_pRecoComparison->Fill(); }

    for (unsigned int i = 0; i < particleVector.size(); ++i) {
      //initialise variables
      m_particle = -1;
      m_pdgcode = 0;
      m_hitsFromSpacePoints = 0;
      m_hitsFromClusters = 0;
      m_hitsFromTrackOrShower = 0;

      const art::Ptr<recob::PFParticle> particle = particleVector.at(i);
      m_particle = particle->Self();
      m_pdgcode = particle->PdgCode();

      PFParticlesToHits::const_iterator pIter = particlesToHits.find(particle);
      PFParticlesToHits::const_iterator pIter2 = particlesToHitsClusters.find(particle);
      if (particlesToHits.end() != pIter) m_hitsFromSpacePoints = pIter->second.size();
      if (particlesToHitsClusters.end() != pIter2) m_hitsFromClusters = pIter2->second.size();

      if (m_pdgcode == 13) {
        PFParticlesToTracks::const_iterator iter = particlesToTracks.find(particle);
        const art::Ptr<recob::Track> track = *(iter->second.begin());

        TracksToHits::const_iterator iter2 = tracksToHits.find(track);
        if (tracksToHits.end() != iter2) {
          const HitVector& hitVector = iter2->second;
          m_hitsFromTrackOrShower = hitVector.size();
        }
      }
      else if (m_pdgcode == 11) {
        PFParticlesToShowers::const_iterator iter = particlesToShowers.find(particle);
        const art::Ptr<recob::Shower> shower = *(iter->second.begin());

        ShowersToHits::const_iterator iter2 = showersToHits.find(shower);
        if (showersToHits.end() != iter2) {
          const HitVector& hitVector = iter2->second;
          m_hitsFromTrackOrShower = hitVector.size();
        }
      }

      if (m_printDebug)
        std::cout << " PFParticle " << m_particle << " (PDG " << m_pdgcode << ") has "
                  << m_hitsFromSpacePoints << " hits from space points and " << m_hitsFromClusters
                  << " hits from clusters, and its recob::Track/Shower has "
                  << m_hitsFromTrackOrShower << " associated hits " << std::endl;

      m_pRecoComparison->Fill();
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void PFParticleHitDumper::FillReco2D(const art::Event& e,
                                       const HitVector& hitVector,
                                       const HitsToPFParticles& hitsToParticles)
  {
    // Initialise variables
    m_particle = -1;
    m_pdgcode = 0;
    m_cstat = 0;
    m_tpc = 0;
    m_plane = 0;
    m_wire = 0;
    m_x = 0.0;
    m_w = 0.0;
    m_q = 0.0;

    // Create dummy entry if there are no 2D hits
    if (hitVector.empty()) { m_pReco2D->Fill(); }

    // Need DetectorProperties service to convert from ticks to X
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);

    // Loop over 2D hits
    for (unsigned int i = 0; i < hitVector.size(); ++i) {
      const art::Ptr<recob::Hit> hit = hitVector.at(i);

      m_particle = -1;
      m_pdgcode = 0;

      HitsToPFParticles::const_iterator pIter = hitsToParticles.find(hit);
      if (hitsToParticles.end() != pIter) {
        const art::Ptr<recob::PFParticle> particle = pIter->second;
        m_particle = particle->Self();
        m_pdgcode = particle->PdgCode();
      }

      const geo::WireID& wireID(hit->WireID());
      m_cstat = wireID.Cryostat;
      m_tpc = wireID.TPC;
      m_plane = wireID.Plane;
      m_wire = wireID.Wire;

      m_q = hit->Integral();
      m_x = detProp.ConvertTicksToX(hit->PeakTime(), wireID.Plane, wireID.TPC, wireID.Cryostat);
      m_w = this->GetUVW(wireID);

      m_pReco2D->Fill();
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void PFParticleHitDumper::FillRecoWires(const art::Event& e, const WireVector& wireVector)
  {

    // Create dummy entry if there are no wires
    if (wireVector.empty()) { m_pRecoWire->Fill(); }

    // Need DetectorProperties service to convert from ticks to X
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);

    // Loop over wires
    int signalCounter(0);

    auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout const>()->Get();
    for (unsigned int i = 0; i < wireVector.size(); ++i) {
      const art::Ptr<recob::Wire> wire = wireVector.at(i);

      const std::vector<float>& signals(wire->Signal());
      const std::vector<geo::WireID> wireIds = wireReadoutGeom.ChannelToWire(wire->Channel());

      if ((signalCounter++) < 10 && m_printDebug)
        std::cout << "    numWires=" << wireVector.size() << " numSignals=" << signals.size()
                  << std::endl;

      double time(0.0);

      m_q = 0.0;

      for (std::vector<float>::const_iterator tIter = signals.begin(), tIterEnd = signals.end();
           tIter != tIterEnd;
           ++tIter) {
        time += 1.0;
        m_q = *tIter;

        if (m_q < 2.0) // seems to remove most noise
          continue;

        for (std::vector<geo::WireID>::const_iterator wIter = wireIds.begin(),
                                                      wIterEnd = wireIds.end();
             wIter != wIterEnd;
             ++wIter) {
          const geo::WireID& wireID = *wIter;
          m_cstat = wireID.Cryostat;
          m_tpc = wireID.TPC;
          m_plane = wireID.Plane;
          m_wire = wireID.Wire;

          m_x = detProp.ConvertTicksToX(time, wireID.Plane, wireID.TPC, wireID.Cryostat);
          m_w = this->GetUVW(wireID);

          m_pRecoWire->Fill();
        }
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  double PFParticleHitDumper::GetUVW(const geo::WireID& wireID) const
  {
    // define UVW as closest distance from (0,0) to wire axis
    auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout>()->Get();

    auto const xyzStart = wireReadoutGeom.Wire(wireID).GetStart();
    const double ay(xyzStart.Y());
    const double az(xyzStart.Z());

    auto const xyzEnd = wireReadoutGeom.Wire(wireID).GetEnd();
    const double by(xyzEnd.Y());
    const double bz(xyzEnd.Z());

    const double ny(by - ay);
    const double nz(bz - az);
    const double N2(ny * ny + nz * nz);

    const double ry(ay - (ay * ny + az * nz) * ny / N2);
    const double rz(az - (ay * ny + az * nz) * nz / N2);
    const double sign((rz > 0.0) ? +1.0 : -1.0);

    return sign * std::sqrt(ry * ry + rz * rz);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  double PFParticleHitDumper::YZtoU(const unsigned int cstat,
                                    const unsigned int tpc,
                                    const double y,
                                    const double z) const
  {
    // TODO: Check that this stills works in DUNE
    auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout>()->Get();
    const double m_theta(wireReadoutGeom.WireAngleToVertical(geo::kU, geo::TPCID{cstat, tpc}));
    return z * std::sin(m_theta) - y * std::cos(m_theta);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  double PFParticleHitDumper::YZtoV(const unsigned int cstat,
                                    const unsigned int tpc,
                                    const double y,
                                    const double z) const
  {
    // TODO; Check that this still works in DUNE
    auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout>()->Get();
    const double m_theta(wireReadoutGeom.WireAngleToVertical(geo::kV, geo::TPCID{cstat, tpc}));
    return z * std::sin(m_theta) - y * std::cos(m_theta);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

} //namespace lar_pandora
