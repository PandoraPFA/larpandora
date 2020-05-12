//############################################################################
//### Name:        ShowerPCAPropergationStartPosition                      ###
//### Author:      Dominic Barker                                          ###
//### Date:        20.09.19                                                ###
//### Description: Get the start position by back propergating the PCA     ###
//###              to the pandora vertex.                                  ###
//############################################################################

#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft Includes
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"

//Root Includes
#include "TVector3.h"

namespace ShowerRecoTools {


  class ShowerPCAPropergationStartPosition: public IShowerTool {

    public:

      ShowerPCAPropergationStartPosition(const fhicl::ParameterSet& pset);

      ~ShowerPCAPropergationStartPosition();

      //Generic Direction Finder
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder
          ) override;

    private:

      //fcl parameters
      art::InputTag fPFParticleModuleLabel;
      std::string   fShowerStartPositionOutputLabel;
      std::string   fShowerCentreInputLabel;
      std::string   fShowerDirectionInputLabel;
      std::string   fShowerStartPositionInputLabel;
  };


  ShowerPCAPropergationStartPosition::ShowerPCAPropergationStartPosition(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel","")),
    fShowerStartPositionOutputLabel(pset.get<std::string>("ShowerStartPositionOutputLabel")),
    fShowerCentreInputLabel(pset.get<std::string>("ShowerCentreInputLabel")),
    fShowerDirectionInputLabel(pset.get<std::string>("ShowerDirectionInputLabel")),
    fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel"))
  {
  }

  ShowerPCAPropergationStartPosition::~ShowerPCAPropergationStartPosition()
  {
  }

  int ShowerPCAPropergationStartPosition::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder){

    TVector3 ShowerCentre        = {-999,-999,-999};

    //Get the start position and direction and center
    if(!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)){
      mf::LogError("ShowerPCAPropergationStartPosition") << "Start position not set, returning "<< std::endl;
      return 1;
    }
    if(!ShowerEleHolder.CheckElement(fShowerDirectionInputLabel)){
      mf::LogError("ShowerPCAPropergationStartPosition") << "Direction not set, returning "<< std::endl;
      return 1;
    }
    if(!ShowerEleHolder.CheckElement(fShowerCentreInputLabel)){

      // Get the assocated pfParicle vertex PFParticles
      art::Handle<std::vector<recob::PFParticle> > pfpHandle;
      if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
        throw cet::exception("ShowerPCADirection") << "Could not get the pandora pf particles. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
        return 1;
      }
      art::FindManyP<recob::SpacePoint>& fmspp = ShowerEleHolder.GetFindManyP<recob::SpacePoint>(
          pfpHandle, Event, fPFParticleModuleLabel);

      if (!fmspp.isValid()){
        throw cet::exception("ShowerPCADirection") << "Trying to get the spacepoint and failed. Something is not configured correctly. Stopping ";
        return 1;
      }

      //Get the spacepoints handle and the hit assoication
      art::Handle<std::vector<recob::SpacePoint> > spHandle;
      if (!Event.getByLabel(fPFParticleModuleLabel, spHandle)){
        throw cet::exception("ShowerPCADirection") << "Could not configure the spacepoint handle. Something is configured incorrectly. Stopping";
        return 1;
      }
      art::FindManyP<recob::Hit>& fmh = ShowerEleHolder.GetFindManyP<recob::Hit>(
          spHandle, Event, fPFParticleModuleLabel);
      if(!fmh.isValid()){
        throw cet::exception("ShowerPCADirection") << "Spacepoint and hit association not valid. Stopping.";
        return 1;
      }

      //Spacepoints
      std::vector<art::Ptr<recob::SpacePoint> > spacePoints_pfp = fmspp.at(pfparticle.key());

      //We cannot progress with no spacepoints.
      if(spacePoints_pfp.size() == 0){return 0;}

      //Get the shower center
      ShowerCentre = IShowerTool::GetLArPandoraShowerAlg().ShowerCentre(spacePoints_pfp,fmh);

    }
    else{
      ShowerEleHolder.GetElement(fShowerCentreInputLabel,ShowerCentre);
    }


    TVector3 ShowerStartPosition = {-999,-999,-999};
    ShowerEleHolder.GetElement(fShowerStartPositionInputLabel,ShowerStartPosition);

    TVector3 ShowerDirection     = {-999,-999,-999};
    ShowerEleHolder.GetElement(fShowerDirectionInputLabel,ShowerDirection);

    //Get the projection
    double projection = ShowerDirection.Dot(ShowerStartPosition-ShowerCentre);

    //Get the position.
    TVector3 ShowerNewStartPosition = projection*ShowerDirection + ShowerCentre;
    TVector3 ShowerNewStartPositionErr = {-999,-999,-999};

    ShowerEleHolder.SetElement(ShowerNewStartPosition,ShowerNewStartPositionErr,fShowerStartPositionOutputLabel);

    return 0;

  }

}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerPCAPropergationStartPosition)

