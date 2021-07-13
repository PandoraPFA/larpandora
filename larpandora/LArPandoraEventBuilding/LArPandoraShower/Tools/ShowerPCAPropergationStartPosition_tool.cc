//############################################################################
//### Name:        ShowerPCAPropergationStartPosition                      ###
//### Author:      Dominic Barker                                          ###
//### Date:        20.09.19                                                ###
//### Description: Get the start position by back propergating the PCA     ###
//###              to the pandora vertex.                                  ###
//############################################################################

//Framework Includes
#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"

namespace ShowerRecoTools {

  class ShowerPCAPropergationStartPosition : public IShowerTool {

  public:
    ShowerPCAPropergationStartPosition(const fhicl::ParameterSet& pset);

    //Generic Direction Finder
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerEleHolder) override;

  private:
    //fcl parameters
    art::InputTag fPFParticleLabel;
    int fVerbose;
    std::string fShowerStartPositionOutputLabel;
    std::string fShowerCentreInputLabel;
    std::string fShowerDirectionInputLabel;
    std::string fShowerStartPositionInputLabel;
  };

  ShowerPCAPropergationStartPosition::ShowerPCAPropergationStartPosition(
    const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel"))
    , fVerbose(pset.get<int>("Verbose"))
    , fShowerStartPositionOutputLabel(pset.get<std::string>("ShowerStartPositionOutputLabel"))
    , fShowerCentreInputLabel(pset.get<std::string>("ShowerCentreInputLabel"))
    , fShowerDirectionInputLabel(pset.get<std::string>("ShowerDirectionInputLabel"))
    , fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel"))
  {}

  int
  ShowerPCAPropergationStartPosition::CalculateElement(
    const art::Ptr<recob::PFParticle>& pfparticle,
    art::Event& Event,
    reco::shower::ShowerElementHolder& ShowerEleHolder)
  {

    TVector3 ShowerCentre = {-999, -999, -999};

    //Get the start position and direction and center
    if (!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)) {
      if (fVerbose)
        mf::LogError("ShowerPCAPropergationStartPosition")
          << "Start position not set, returning " << std::endl;
      return 1;
    }
    if (!ShowerEleHolder.CheckElement(fShowerDirectionInputLabel)) {
      if (fVerbose)
        mf::LogError("ShowerPCAPropergationStartPosition")
          << "Direction not set, returning " << std::endl;
      return 1;
    }
    if (!ShowerEleHolder.CheckElement(fShowerCentreInputLabel)) {

      auto const clockData =
        art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(Event);
      auto const detProp =
        art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(Event, clockData);

      // Get the assocated pfParicle vertex PFParticles
      auto const pfpHandle = Event.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);

      const art::FindManyP<recob::SpacePoint>& fmspp =
        ShowerEleHolder.GetFindManyP<recob::SpacePoint>(pfpHandle, Event, fPFParticleLabel);

      //Get the spacepoints handle and the hit assoication
      auto const spHandle = Event.getValidHandle<std::vector<recob::SpacePoint>>(fPFParticleLabel);

      const art::FindManyP<recob::Hit>& fmh =
        ShowerEleHolder.GetFindManyP<recob::Hit>(spHandle, Event, fPFParticleLabel);

      //Spacepoints
      std::vector<art::Ptr<recob::SpacePoint>> spacePoints_pfp = fmspp.at(pfparticle.key());

      //We cannot progress with no spacepoints.
      if (spacePoints_pfp.empty()) return 1;

      //Get the shower center
      ShowerCentre = IShowerTool::GetLArPandoraShowerAlg().ShowerCentre(
        clockData, detProp, spacePoints_pfp, fmh);
    }
    else {
      ShowerEleHolder.GetElement(fShowerCentreInputLabel, ShowerCentre);
    }

    TVector3 ShowerStartPosition = {-999, -999, -999};
    ShowerEleHolder.GetElement(fShowerStartPositionInputLabel, ShowerStartPosition);

    TVector3 ShowerDirection = {-999, -999, -999};
    ShowerEleHolder.GetElement(fShowerDirectionInputLabel, ShowerDirection);

    //Get the projection
    double projection = ShowerDirection.Dot(ShowerStartPosition - ShowerCentre);

    //Get the position.
    TVector3 ShowerNewStartPosition = projection * ShowerDirection + ShowerCentre;
    TVector3 ShowerNewStartPositionErr = {-999, -999, -999};

    ShowerEleHolder.SetElement(
      ShowerNewStartPosition, ShowerNewStartPositionErr, fShowerStartPositionOutputLabel);

    return 0;
  }

}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerPCAPropergationStartPosition)
