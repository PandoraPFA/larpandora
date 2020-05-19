//############################################################################
//### Name:        ShowerTrackStartPosition                                ###
//### Author:      You                                                     ###
//### Date:        13.05.19                                                ###
//### Description: Generic form of the shower tools                        ###
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
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"

namespace ShowerRecoTools {


  class ShowerTrackStartPosition: public IShowerTool {

    public:

      ShowerTrackStartPosition(const fhicl::ParameterSet& pset);

      ~ShowerTrackStartPosition();

      //Generic Direction Finder
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder
          ) override;

    private:

      int         fVerbose;
      std::string fInitialTrackInputLabel;
      std::string fShowerStartPositionOutputLabel;

  };


  ShowerTrackStartPosition::ShowerTrackStartPosition(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fVerbose(pset.get<int>("Verbose")),
    fInitialTrackInputLabel(pset.get<std::string>("InitialTrackInputLabel")),
    fShowerStartPositionOutputLabel(pset.get<std::string>("ShowerStartPositionOutputLabel"))
  {
  }

  ShowerTrackStartPosition::~ShowerTrackStartPosition()
  {
  }

  int ShowerTrackStartPosition::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder){

    //Check the Track has been defined
    if(!ShowerEleHolder.CheckElement("InitialTrack")){
      if (fVerbose)
        mf::LogError("ShowerSmartTrackTrajectoryPointDirection")
          << "Initial track not set"<< std::endl;
      return 1;
    }
    recob::Track InitialTrack;
    ShowerEleHolder.GetElement(fInitialTrackInputLabel,InitialTrack);


    //Set the shower start position as the
    TVector3 StartPositionErr = {-999,-999,-999};

    geo::Point_t  TrajPosition_vec   = InitialTrack.LocationAtPoint(0);
    TVector3 TrajPosition = {TrajPosition_vec.X(), TrajPosition_vec.Y(),TrajPosition_vec.Z()};
    ShowerEleHolder.SetElement(TrajPosition,StartPositionErr,fShowerStartPositionOutputLabel);

    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerTrackStartPosition)

