//############################################################################
//### Name:        ShowerLengthPercentile                                  ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Simple code to calculate the lenght such that a given % ###
//###              of the hits are within the length.                      ###
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

namespace ShowerRecoTools {


  class ShowerLengthPercentile: public IShowerTool {

    public:

      ShowerLengthPercentile(const fhicl::ParameterSet& pset);

      //Generic Direction Finder
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder
          ) override;

    private:

      float fPercentile;

      art::InputTag fPFParticleLabel;
      int           fVerbose;
      std::string   fShowerStartPositionInputLabel;
      std::string   fShowerDirectionInputLabel;
      std::string   fShowerLengthOutputLabel;
      std::string   fShowerOpeningAngleOutputLabel;
  };


  ShowerLengthPercentile::ShowerLengthPercentile(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fPercentile(pset.get<float>("Percentile")),
    fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel")),
    fVerbose(pset.get<int>("Verbose")),
    fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel")),
    fShowerDirectionInputLabel(pset.get<std::string>("ShowerDirectionInputLabel")),
    fShowerLengthOutputLabel(pset.get<std::string>("ShowerLengthOutputLabel")),
    fShowerOpeningAngleOutputLabel(pset.get<std::string>("ShowerOpeningAngleOutputLabel"))
  {
  }

  int ShowerLengthPercentile::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder){

    //Get the start position
    if(!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)){
      if (fVerbose)
        mf::LogError("ShowerSlidingStandardCalodEdx") << "Start position not set, returning "<< std::endl;
      return 1;
    }
    //Only consider hits in the same tpcs as the vertex.
    TVector3 ShowerStartPosition = {-999,-999,-999};
    ShowerEleHolder.GetElement(fShowerStartPositionInputLabel,ShowerStartPosition);

    // Get the assocated pfParicle Handle
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleLabel, pfpHandle)){
      throw cet::exception("ShowerResidualTrackHitFinder") << "Could not get the pandora pf particles. Something is not cofingured correctly Please give the correct pandoa module label. Stopping";
      return 1;
    }

    // Get the spacepoint - PFParticle assn
    art::FindManyP<recob::SpacePoint>& fmspp = ShowerEleHolder.GetFindManyP<recob::SpacePoint>(
        pfpHandle, Event, fPFParticleLabel);
    if (!fmspp.isValid()){
      throw cet::exception("ShowerResidualTrackHitFinder") << "Trying to get the spacepoint and failed. Something is not configured correctly. Stopping ";
      return 1;
    }

    // Get the spacepoints
    art::Handle<std::vector<recob::SpacePoint> > spHandle;
    if (!Event.getByLabel(fPFParticleLabel, spHandle)){
      throw cet::exception("ShowerResidualTrackHitFinder") << "Could not configure the spacepoint handle. Something is configured incorrectly. Stopping";
      return 1;
    }

    // Get the SpacePoints
    std::vector<art::Ptr<recob::SpacePoint> > spacePoints = fmspp.at(pfparticle.key());
    if (!spacePoints.size()){
      if (fVerbose)
        mf::LogError("ShowerLengthPercentile") << "No Spacepoints, returning" <<std::endl;
      return 1;
    }


    if(!ShowerEleHolder.CheckElement(fShowerDirectionInputLabel)){
      if (fVerbose)
        mf::LogError("ShowerResidualTrackHitFinder") << "Direction not set, returning "<< std::endl;
      return 1;
    }

    TVector3 ShowerDirection     = {-999,-999,-999};
    ShowerEleHolder.GetElement(fShowerDirectionInputLabel,ShowerDirection);

    //Order the spacepoints
    IShowerTool::GetLArPandoraShowerAlg().OrderShowerSpacePoints(spacePoints,ShowerStartPosition,ShowerDirection);

    //Find the length as the value that contains % of the hits
    int lengthIter = fPercentile*spacePoints.size();

    //Find the length
    double ShowerLength = IShowerTool::GetLArPandoraShowerAlg().SpacePointProjection(
        spacePoints[lengthIter], ShowerStartPosition, ShowerDirection);
    double ShowerMaxProjection = IShowerTool::GetLArPandoraShowerAlg().SpacePointProjection(
        spacePoints[spacePoints.size() -1], ShowerStartPosition, ShowerDirection);

    double ShowerLengthError = ShowerMaxProjection - ShowerLength;

    //Order the spacepoints in perpendicular
    IShowerTool::GetLArPandoraShowerAlg().OrderShowerSpacePointsPerpendicular(spacePoints,ShowerStartPosition,ShowerDirection);

    //Find the length as the value that contains % of the hits
    int perpIter = fPercentile*spacePoints.size();

    //Find the width of the shower
    double ShowerWidth = IShowerTool::GetLArPandoraShowerAlg().SpacePointPerpendicular(
        spacePoints[perpIter], ShowerStartPosition, ShowerDirection);
    // double ShowerMaxWidth = IShowerTool::GetLArPandoraShowerAlg().SpacePointPerpendicular(
    //     spacePoints[spacePoints.size() -1], ShowerStartPosition, ShowerDirection);

    double ShowerAngle = atan(ShowerWidth/ShowerLength);
    double ShowerAngleError = -9999; //TODO: Do properly

    // Fill the shower element holder
    ShowerEleHolder.SetElement(ShowerLength, ShowerLengthError, fShowerLengthOutputLabel);
    ShowerEleHolder.SetElement(ShowerAngle, ShowerAngleError, fShowerOpeningAngleOutputLabel);

    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerLengthPercentile)
