//############################################################################
//### Name:        ShowerTrackHitDirection                                 ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the shower direction using the         ###
//###              initial track the average direction of the initial hits ###
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
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"

//C++ Includes
#include <iostream>

//Root Includes
#include "TVector3.h"
#include "TMath.h"

namespace ShowerRecoTools {


  class ShowerTrackHitDirection:IShowerTool {

    public:

      ShowerTrackHitDirection(const fhicl::ParameterSet& pset);

      ~ShowerTrackHitDirection();

      //Calculate the shower direction from the initial track hits.
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder
          ) override;

    private:

      //fcl
      int  fVerbose;
      bool fUsePandoraVertex; //Direction from point defined as (Position of Hit - Vertex)
      //rather than (Position of Hit - Track Start Point)
      art::InputTag fHitModuleLabel;
      art::InputTag fPFParticleLabel;

      std::string fInitialTrackHitsInputLabel;
      std::string fShowerStartPositionInputLabel;
      std::string fInitialTrackInputLabel;
      std::string fShowerDirectionOutputLabel;
  };


  ShowerTrackHitDirection::ShowerTrackHitDirection(const fhicl::ParameterSet& pset)
    :  IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fVerbose(pset.get<int>("Verbose")),
    fUsePandoraVertex(pset.get<bool>         ("UsePandoraVertex")),
    fHitModuleLabel(pset.get<art::InputTag>("HitModuleLabel")),
    fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel")),
    fInitialTrackHitsInputLabel(pset.get<std::string>("InitialTrackHitsInputLabel")),
    fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel")),
    fInitialTrackInputLabel(pset.get<std::string>("InitialTrackInputLabel")),
    fShowerDirectionOutputLabel(pset.get<std::string>("ShowerDirectionOutputLabel"))
  {
  }

  ShowerTrackHitDirection::~ShowerTrackHitDirection()
  {
  }

  int ShowerTrackHitDirection::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event,
      reco::shower::ShowerElementHolder& ShowerEleHolder){

    //Check the Track Hits has been defined
    if(!ShowerEleHolder.CheckElement(fInitialTrackHitsInputLabel)){
      if (fVerbose)
        mf::LogError("ShowerTrackHitDirection") << "Initial track hits not set"<< std::endl;
      return 0;
    }

    //Check the start position is set.
    if(fUsePandoraVertex && !ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)){
      if (fVerbose)
        mf::LogError("ShowerTrackHitDirection") << "Start position not set, returning "<< std::endl;
      return 0;
    }

    //Get the start poistion
    TVector3 StartPosition = {-999,-999,-99};
    if(fUsePandoraVertex){
      ShowerEleHolder.GetElement(fShowerStartPositionInputLabel,StartPosition);
    }
    else{
      //Check the Tracks has been defined
      if(!ShowerEleHolder.CheckElement(fInitialTrackInputLabel)){
        if (fVerbose)
          mf::LogError("ShowerTrackHitDirection") << "Initial track not set"<< std::endl;
        return 0;
      }
      recob::Track InitialTrack;
      ShowerEleHolder.GetElement(fInitialTrackInputLabel,InitialTrack);
      geo::Point_t Start_point = InitialTrack.Start();
      StartPosition = {Start_point.X(),Start_point.Y(),Start_point.Z()};
    }


    //Get the spacepoints associated to hits
    art::Handle<std::vector<recob::Hit> > hitHandle;
    if (!Event.getByLabel(fHitModuleLabel, hitHandle)){
      throw cet::exception("ShowerTrackHitDirection") << "Could not configure the spacepoint handle. Something is configured incorrectly. Stopping";
      return 1;
    }

    //Get the spacepoint handle. We need to do this in 3D.
    art::FindManyP<recob::SpacePoint>& fmsp = ShowerEleHolder.GetFindManyP<recob::SpacePoint>(
        hitHandle, Event, fPFParticleLabel);
    if(!fmsp.isValid()){
      throw cet::exception("ShowerTrackHitDirection") << "Spacepoint and hit association not valid. Stopping.";
      return 1;
    }

    //Get the initial track hits.
    std::vector<art::Ptr<recob::Hit> > InitialTrackHits;
    ShowerEleHolder.GetElement(fInitialTrackHitsInputLabel,InitialTrackHits);


    //Calculate the mean direction and the the standard deviation
    float sumX=0, sumX2=0;
    float sumY=0, sumY2=0;
    float sumZ=0, sumZ2=0;

    //Get the spacepoints associated to the track hit
    std::vector<art::Ptr<recob::SpacePoint > > intitaltrack_sp;
    for(auto const hit: InitialTrackHits){
      std::vector<art::Ptr<recob::SpacePoint > > sps = fmsp.at(hit.key());
      for(auto const sp: sps){
        intitaltrack_sp.push_back(sp);

        //Get the direction relative to the start positon
        TVector3 pos = IShowerTool::GetLArPandoraShowerAlg().SpacePointPosition(sp) - StartPosition;
        if(pos.Mag() == 0){continue;}

        sumX = pos.X(); sumX2 += pos.X()*pos.X();
        sumY = pos.Y(); sumY2 += pos.Y()*pos.Y();
        sumZ = pos.Z(); sumZ2 += pos.Z()*pos.Z();
      }
    }

    float NumSps = intitaltrack_sp.size();
    TVector3 Mean = {sumX/NumSps,sumY/NumSps,sumZ/NumSps};
    Mean = Mean.Unit();

    float RMSX = 999;
    float RMSY = 999;
    float RMSZ = 999;
    if(sumX2/NumSps - ((sumX/NumSps)*((sumX/NumSps))) > 0){
      RMSX = TMath::Sqrt(sumX2/NumSps - ((sumX/NumSps)*((sumX/NumSps))));
    }
    if(sumY2/NumSps - ((sumY/NumSps)*((sumY/NumSps))) > 0){
      RMSY = TMath::Sqrt(sumY2/NumSps - ((sumY/NumSps)*((sumY/NumSps))));
    }
    if(sumZ2/NumSps - ((sumZ/NumSps)*((sumZ/NumSps))) > 0){
      RMSZ = TMath::Sqrt(sumZ2/NumSps - ((sumZ/NumSps)*((sumZ/NumSps))));
    }


    //Loop over the spacepoints and remove ones the relative direction is not within one sigma.
    TVector3 Direction_Mean = {0,0,0};
    int N = 0;
    for(auto const sp: intitaltrack_sp){
      TVector3 Direction = IShowerTool::GetLArPandoraShowerAlg().SpacePointPosition(sp) - StartPosition;
      if((TMath::Abs((Direction-Mean).X()) < 1*RMSX) &&
          (TMath::Abs((Direction-Mean).Y())< 1*RMSY) &&
          (TMath::Abs((Direction-Mean).Z()) < 1*RMSZ)){
        if(Direction.Mag() == 0){continue;}
        ++N;
        Direction_Mean += Direction;
      }
    }

    if(N>0){
      //Take the mean value
      TVector3 Direction = Direction_Mean.Unit();
      ShowerEleHolder.SetElement(Direction,fShowerDirectionOutputLabel);
    }
    else{
      if (fVerbose)
        mf::LogError("ShowerTrackHitDirection") << "None of the points are within 1 sigma"<< std::endl;
      return 1;
    }
    return 0;
  }
}


DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerTrackHitDirection)

