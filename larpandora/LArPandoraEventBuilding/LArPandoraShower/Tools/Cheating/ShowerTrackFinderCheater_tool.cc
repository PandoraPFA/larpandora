//############################################################################
//### Name:        ShowerTrackFinderCheater                                ###
//### Author:      Ed Tyley                                                ###
//### Date:        16.07.19                                                ###
//### Description: Cheating tool using truth for shower direction          ###
//############################################################################

//Framework Includes
#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "lardataobj/RecoBase/Cluster.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Algs/LArPandoraShowerCheatingAlg.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"

namespace ShowerRecoTools {

  class ShowerTrackFinderCheater : IShowerTool {

  public:
    ShowerTrackFinderCheater(const fhicl::ParameterSet& pset);

    //Generic Direction Finder
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerEleHolder) override;

  private:
    //Algorithm functions
    shower::LArPandoraShowerCheatingAlg fLArPandoraShowerCheatingAlg;

    //fcl
    const bool fDebugEVD;
    const art::InputTag fPFParticleLabel;
    const art::InputTag fHitModuleLabel;

    const std::string fTrueParticleIntputLabel;
    const std::string fShowerStartPositionInputTag;
    const std::string fShowerDirectionInputTag;
    const std::string fInitialTrackHitsOutputLabel;
    const std::string fInitialTrackSpacePointsOutputLabel;
  };

  ShowerTrackFinderCheater::ShowerTrackFinderCheater(const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fLArPandoraShowerCheatingAlg(pset.get<fhicl::ParameterSet>("LArPandoraShowerCheatingAlg"))
    , fDebugEVD(pset.get<bool>("DebugEVD"))
    , fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel"))
    , fHitModuleLabel(pset.get<art::InputTag>("HitModuleLabel"))
    , fTrueParticleIntputLabel(pset.get<std::string>("TrueParticleIntputLabel"))
    , fShowerStartPositionInputTag(pset.get<std::string>("ShowerStartPositionInputTag"))
    , fShowerDirectionInputTag(pset.get<std::string>("ShowerDirectionInputTag"))
    , fInitialTrackHitsOutputLabel(pset.get<std::string>("InitialTrackHitsOutputLabel"))
    , fInitialTrackSpacePointsOutputLabel(
        pset.get<std::string>("InitialTrackSpacePointsOutputLabel"))
  {}

  int
  ShowerTrackFinderCheater::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                                             art::Event& Event,
                                             reco::shower::ShowerElementHolder& ShowerEleHolder)
  {

    const simb::MCParticle* trueParticle;

    //Get the hits from the shower:
    auto const pfpHandle = Event.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
    auto const clockData =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(Event);

    //Get the clusters
    auto const clusHandle = Event.getValidHandle<std::vector<recob::Cluster>>(fPFParticleLabel);
    art::FindManyP<recob::Cluster> fmc(pfpHandle, Event, fPFParticleLabel);
    std::vector<art::Ptr<recob::Cluster>> clusters = fmc.at(pfparticle.key());

    //Get the hit association
    art::FindManyP<recob::Hit> fmhc(clusHandle, Event, fPFParticleLabel);

    std::vector<art::Ptr<recob::Hit>> showerHits;
    for (auto const& cluster : clusters) {

      //Get the hits
      std::vector<art::Ptr<recob::Hit>> hits = fmhc.at(cluster.key());
      showerHits.insert(showerHits.end(), hits.begin(), hits.end());
    }

    if (ShowerEleHolder.CheckElement(fTrueParticleIntputLabel)) {
      ShowerEleHolder.GetElement(fTrueParticleIntputLabel, trueParticle);
    }
    else {

      //Could store these in the shower element holder and just calculate once?
      std::map<int, const simb::MCParticle*> trueParticles =
        fLArPandoraShowerCheatingAlg.GetTrueParticleMap();
      std::map<int, std::vector<int>> showersMothers =
        fLArPandoraShowerCheatingAlg.GetTrueChain(trueParticles);

      //Get the true particle from the shower
      std::pair<int, double> ShowerTrackInfo =
        fLArPandoraShowerCheatingAlg.TrueParticleIDFromTrueChain(
          clockData, showersMothers, showerHits, 2);

      if (ShowerTrackInfo.first == -99999) {
        mf::LogError("ShowerStartPosition") << "True Shower Not Found";
        return 1;
      }
      trueParticle = trueParticles[ShowerTrackInfo.first];
      ShowerEleHolder.SetElement(trueParticle, fTrueParticleIntputLabel);
    }

    if (!trueParticle) {
      mf::LogError("ShowerDirectionCheater") << "True shower not found, returning";
      return 1;
    }

    //This is all based on the shower vertex being known. If it is not lets not do the track
    if (!ShowerEleHolder.CheckElement(fShowerStartPositionInputTag)) {
      mf::LogError("ShowerTrackFinderCheater") << "Start position not set, returning " << std::endl;
      return 1;
    }
    if (!ShowerEleHolder.CheckElement(fShowerDirectionInputTag)) {
      mf::LogError("ShowerTrackFinderCheater") << "Direction not set, returning " << std::endl;
      return 1;
    }

    TVector3 ShowerStartPosition = {-999, -999, -999};
    ShowerEleHolder.GetElement(fShowerStartPositionInputTag, ShowerStartPosition);

    TVector3 ShowerDirection = {-999, -999, -999};
    ShowerEleHolder.GetElement(fShowerDirectionInputTag, ShowerDirection);

    auto const hitHandle = Event.getValidHandle<std::vector<recob::Hit>>(fHitModuleLabel);

    // Get the hits associated with the space points
    art::FindManyP<recob::SpacePoint> fmsph(hitHandle, Event, fPFParticleLabel);
    if (!fmsph.isValid()) {
      throw cet::exception("ShowerTrackFinderCheater")
        << "Spacepoint and hit association not valid. Stopping.";
    }

    std::vector<int> trueParticleIdVec;

    // If we have an electron, take the hits from the primary only
    // This will also cover the cases Pandora misclassify a track as a shower
    if (trueParticle->PdgCode() != 22) { trueParticleIdVec.push_back(trueParticle->TrackId()); }
    else {
      // To check if we are rolling up showers, check the number of daughters the photon has
      const int nDaughters = trueParticle->NumberDaughters();
      if (nDaughters == 0) {
        // If we roll up showers, we have no choice but to take all of the hits from the photon
        trueParticleIdVec.push_back(-trueParticle->TrackId());
      }
      else {
        // If we do not roll up the showers, take all of the primary daughters
        for (int i = 0; i < nDaughters; i++) {
          trueParticleIdVec.push_back(trueParticle->Daughter(i));
        }
      }
    }

    std::vector<art::Ptr<recob::Hit>> trackHits;
    std::vector<art::Ptr<recob::SpacePoint>> trackSpacePoints;

    //Get the hits from the true particle
    for (auto hit : showerHits) {
      int trueHitId = fLArPandoraShowerCheatingAlg.TrueParticleID(clockData, hit);
      if (std::find(trueParticleIdVec.cbegin(), trueParticleIdVec.cend(), trueHitId) !=
          trueParticleIdVec.cend()) {
        trackHits.push_back(hit);
        std::vector<art::Ptr<recob::SpacePoint>> sps = fmsph.at(hit.key());
        if (sps.size() == 1) { trackSpacePoints.push_back(sps.front()); }
      }
    }

    if (trackHits.empty() || trackSpacePoints.empty())
      mf::LogWarning("ShowerTrackFinderCheater")
        << "Creating intial track with " << trackHits.size() << " hits and "
        << trackSpacePoints.size() << " spacepoints" << std::endl;

    ShowerEleHolder.SetElement(trackHits, fInitialTrackHitsOutputLabel);
    ShowerEleHolder.SetElement(trackSpacePoints, fInitialTrackSpacePointsOutputLabel);

    if (fDebugEVD) {
      fLArPandoraShowerCheatingAlg.CheatDebugEVD(
        clockData, trueParticle, Event, ShowerEleHolder, pfparticle);
    }

    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerTrackFinderCheater)
