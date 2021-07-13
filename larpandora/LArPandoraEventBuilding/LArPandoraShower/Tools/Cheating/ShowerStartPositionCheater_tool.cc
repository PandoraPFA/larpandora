//############################################################################
//### Name:        ShowerStartPositionCheater                              ###
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

  class ShowerStartPositionCheater : IShowerTool {

  public:
    ShowerStartPositionCheater(const fhicl::ParameterSet& pset);

    //Calculate Cheating Start Position
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerEleHolder) override;

  private:
    //Algorithm functions
    shower::LArPandoraShowerCheatingAlg fLArPandoraShowerCheatingAlg;

    //FCL
    const art::InputTag fPFParticleLabel;
    const art::InputTag fHitModuleLabel;

    const std::string fShowerStartPositionOutputLabel;
    const std::string fTrueParticleOutputLabel;
  };

  ShowerStartPositionCheater::ShowerStartPositionCheater(const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fLArPandoraShowerCheatingAlg(pset.get<fhicl::ParameterSet>("LArPandoraShowerCheatingAlg"))
    , fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel"))
    , fHitModuleLabel(pset.get<art::InputTag>("HitModuleLabel"))
    , fShowerStartPositionOutputLabel(pset.get<std::string>("ShowerStartPositionOutputLabel"))
    , fTrueParticleOutputLabel(pset.get<std::string>("TrueParticleOutputLabel"))
  {}

  int
  ShowerStartPositionCheater::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                                               art::Event& Event,
                                               reco::shower::ShowerElementHolder& ShowerEleHolder)
  {

    //Could store these in the shower element holder and just calculate once?
    std::map<int, const simb::MCParticle*> trueParticles =
      fLArPandoraShowerCheatingAlg.GetTrueParticleMap();
    std::map<int, std::vector<int>> showersMothers =
      fLArPandoraShowerCheatingAlg.GetTrueChain(trueParticles);

    //Get the hits from the shower:
    auto const pfpHandle = Event.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);

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

    //Get the true particle from the shower
    auto const clockData =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(Event);

    std::pair<int, double> ShowerTrackInfo =
      fLArPandoraShowerCheatingAlg.TrueParticleIDFromTrueChain(
        clockData, showersMothers, showerHits, 2);

    if (ShowerTrackInfo.first == -99999) {
      mf::LogError("ShowerStartPositionCheater") << "True Shower Not Found";
      return 1;
    }

    const simb::MCParticle* trueParticle = trueParticles[ShowerTrackInfo.first];
    if (!trueParticle) {
      mf::LogError("ShowerDirectionCheater") << "True shower not found, returning";
      return 1;
    }

    ShowerEleHolder.SetElement(trueParticle, fTrueParticleOutputLabel);

    TVector3 trueStartPos = {-999, -999, -999};
    // If the true particle is a photon, we need to be smarter.
    // Select the first traj point where tne photon loses energy
    if (abs(trueParticle->PdgCode()) == 22) {
      double initialEnergy = trueParticle->E();
      unsigned int TrajPoints = trueParticle->NumberTrajectoryPoints();
      for (unsigned int trajPoint = 0; trajPoint < TrajPoints; trajPoint++) {
        if (trueParticle->E(trajPoint) < initialEnergy) {
          trueStartPos = trueParticle->Position(trajPoint).Vect();
          break;
        }
      }
    }
    else {
      trueStartPos = trueParticle->Position().Vect();
    }

    TVector3 trueStartPosErr = {-999, -999, -999};
    ShowerEleHolder.SetElement(trueStartPos, trueStartPosErr, fShowerStartPositionOutputLabel);

    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerStartPositionCheater)
