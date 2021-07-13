//############################################################################
//### Name:        ShowerDirectionCheater                                  ###
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

  class ShowerDirectionCheater : IShowerTool {

  public:
    ShowerDirectionCheater(const fhicl::ParameterSet& pset);

    //Generic Direction Finder
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerEleHolder) override;

  private:
    //Algorithm functions
    shower::LArPandoraShowerCheatingAlg fLArPandoraShowerCheatingAlg;

    //Services
    art::ServiceHandle<art::TFileService> tfs;

    //fcl
    const art::InputTag fPFParticleLabel;
    const unsigned int
      fNSegments;        //Number of segement to split the shower into the perforam the RMSFlip.
    const bool fRMSFlip; //Flip the direction by considering the rms.
    const bool
      fVertexFlip; //Flip the direction by considering the vertex position relative to the center position.

    //TTree Branch variables
    TTree* Tree;
    float vertexDotProduct;
    float rmsGradient;

    const std::string fShowerStartPositionInputLabel;
    const std::string fTrueParticleInputLabel;
    const std::string fShowerDirectionOutputLabel;
  };

  ShowerDirectionCheater::ShowerDirectionCheater(const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fLArPandoraShowerCheatingAlg(pset.get<fhicl::ParameterSet>("LArPandoraShowerCheatingAlg"))
    , fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel"))
    , fNSegments(pset.get<unsigned int>("NSegments"))
    , fRMSFlip(pset.get<bool>("RMSFlip"))
    , fVertexFlip(pset.get<bool>("VertexFlip"))
    , fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel"))
    , fTrueParticleInputLabel(pset.get<std::string>("TrueParticleInputLabel"))
    , fShowerDirectionOutputLabel(pset.get<std::string>("ShowerDirectionOutputLabel"))
  {
    if (fVertexFlip || fRMSFlip) {
      Tree = tfs->make<TTree>("DebugTreeDirCheater", "DebugTree from shower direction cheater");
      if (fVertexFlip) Tree->Branch("vertexDotProduct", &vertexDotProduct);
      if (fRMSFlip) Tree->Branch("rmsGradient", &rmsGradient);
    }
  }

  int
  ShowerDirectionCheater::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                                           art::Event& Event,
                                           reco::shower::ShowerElementHolder& ShowerEleHolder)
  {

    const simb::MCParticle* trueParticle;

    //Get the hits from the shower:
    auto const pfpHandle = Event.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);

    auto const clockData =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(Event);
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(Event, clockData);

    if (ShowerEleHolder.CheckElement(fTrueParticleInputLabel)) {
      ShowerEleHolder.GetElement(fTrueParticleInputLabel, trueParticle);
    }
    else {

      //Could store these in the shower element holder and just calculate once?
      std::map<int, const simb::MCParticle*> trueParticles =
        fLArPandoraShowerCheatingAlg.GetTrueParticleMap();
      std::map<int, std::vector<int>> showersMothers =
        fLArPandoraShowerCheatingAlg.GetTrueChain(trueParticles);

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
      std::pair<int, double> ShowerTrackInfo =
        fLArPandoraShowerCheatingAlg.TrueParticleIDFromTrueChain(
          clockData, showersMothers, showerHits, 2);

      if (ShowerTrackInfo.first == -99999) {
        mf::LogError("ShowerDirectionCheater") << "True shower not found, returning";
        return 1;
      }
      trueParticle = trueParticles[ShowerTrackInfo.first];
      ShowerEleHolder.SetElement(trueParticle, fTrueParticleInputLabel);
    }

    if (!trueParticle) {
      mf::LogError("ShowerDirectionCheater") << "True shower not found, returning";
      return 1;
    }

    TVector3 trueDir = TVector3{trueParticle->Px(), trueParticle->Py(), trueParticle->Pz()}.Unit();

    TVector3 trueDirErr = {-999, -999, -999};
    ShowerEleHolder.SetElement(trueDir, trueDirErr, fShowerDirectionOutputLabel);

    if (fRMSFlip || fVertexFlip) {

      // Reset the tree values to defaults
      rmsGradient = std::numeric_limits<float>::lowest();
      vertexDotProduct = std::numeric_limits<float>::lowest();

      //Get the SpacePoints and hits
      art::FindManyP<recob::SpacePoint> fmspp(pfpHandle, Event, fPFParticleLabel);

      if (!fmspp.isValid()) {
        throw cet::exception("ShowerDirectionCheater")
          << "Trying to get the spacepoint and failed. Something is not configured correctly. "
             "Stopping ";
      }

      auto const spHandle = Event.getValidHandle<std::vector<recob::SpacePoint>>(fPFParticleLabel);
      art::FindManyP<recob::Hit> fmh(spHandle, Event, fPFParticleLabel);
      if (!fmh.isValid()) {
        throw cet::exception("ShowerDirectionCheater")
          << "Spacepoint and hit association not valid. Stopping.";
      }
      std::vector<art::Ptr<recob::SpacePoint>> spacePoints = fmspp.at(pfparticle.key());

      if (spacePoints.size() < 3) {
        mf::LogWarning("ShowerDirectionCheater")
          << spacePoints.size() << " spacepoints in shower, not calculating direction" << std::endl;
        return 1;
      }

      //Get Shower Centre
      float TotalCharge;

      const TVector3 ShowerCentre = IShowerTool::GetLArPandoraShowerAlg().ShowerCentre(
        clockData, detProp, spacePoints, fmh, TotalCharge);

      //Check if we are pointing the correct direction or not, First try the start position
      if (ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel) && fVertexFlip) {

        //Get the General direction as the vector between the start position and the centre
        TVector3 StartPositionVec = {-999, -999, -999};
        ShowerEleHolder.GetElement(fShowerStartPositionInputLabel, StartPositionVec);

        TVector3 GeneralDir = (ShowerCentre - StartPositionVec).Unit();

        //Dot product
        vertexDotProduct = trueDir.Dot(GeneralDir);
      }

      if (fRMSFlip) {
        //Otherwise Check against the RMS of the shower. Method adapated from EMShower Thanks Mike.
        rmsGradient = IShowerTool::GetLArPandoraShowerAlg().RMSShowerGradient(
          spacePoints, ShowerCentre, trueDir, fNSegments);
      }
      Tree->Fill();
    }

    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerDirectionCheater)
