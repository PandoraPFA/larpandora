//############################################################################
//### Name:        ShowerLengthPercentile                                  ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Simple code to calculate the lenght such that a given % ###
//###              of the hits are within the length.                      ###
//############################################################################

//Framework Includes
#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"

namespace ShowerRecoTools {

  class ShowerLengthPercentile : public IShowerTool {

  public:
    ShowerLengthPercentile(const fhicl::ParameterSet& pset);

    //Generic Direction Finder
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerEleHolder) override;

  private:
    float fPercentile;

    art::InputTag fPFParticleLabel;
    int fVerbose;
    std::string fShowerStartPositionInputLabel;
    std::string fShowerDirectionInputLabel;
    std::string fShowerLengthOutputLabel;
    std::string fShowerOpeningAngleOutputLabel;
  };

  ShowerLengthPercentile::ShowerLengthPercentile(const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fPercentile(pset.get<float>("Percentile"))
    , fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel"))
    , fVerbose(pset.get<int>("Verbose"))
    , fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel"))
    , fShowerDirectionInputLabel(pset.get<std::string>("ShowerDirectionInputLabel"))
    , fShowerLengthOutputLabel(pset.get<std::string>("ShowerLengthOutputLabel"))
    , fShowerOpeningAngleOutputLabel(pset.get<std::string>("ShowerOpeningAngleOutputLabel"))
  {}

  int
  ShowerLengthPercentile::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                                           art::Event& Event,
                                           reco::shower::ShowerElementHolder& ShowerEleHolder)
  {

    //Get the start position
    if (!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)) {
      if (fVerbose)
        mf::LogError("ShowerLengthPercentile") << "Start position not set, returning " << std::endl;
      return 1;
    }
    //Only consider hits in the same tpcs as the vertex.
    TVector3 ShowerStartPosition = {-999, -999, -999};
    ShowerEleHolder.GetElement(fShowerStartPositionInputLabel, ShowerStartPosition);

    // Get the assocated pfParicle Handle
    auto const pfpHandle = Event.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);

    // Get the spacepoint - PFParticle assn
    const art::FindManyP<recob::SpacePoint>& fmspp =
      ShowerEleHolder.GetFindManyP<recob::SpacePoint>(pfpHandle, Event, fPFParticleLabel);

    // Get the SpacePoints
    std::vector<art::Ptr<recob::SpacePoint>> spacePoints = fmspp.at(pfparticle.key());
    if (spacePoints.empty()) {
      if (fVerbose)
        mf::LogError("ShowerLengthPercentile") << "No Spacepoints, returning" << std::endl;
      return 1;
    }

    if (!ShowerEleHolder.CheckElement(fShowerDirectionInputLabel)) {
      if (fVerbose)
        mf::LogError("ShowerLengthPercentile") << "Direction not set, returning " << std::endl;
      return 1;
    }

    TVector3 ShowerDirection = {-999, -999, -999};
    ShowerEleHolder.GetElement(fShowerDirectionInputLabel, ShowerDirection);

    //Order the spacepoints
    IShowerTool::GetLArPandoraShowerAlg().OrderShowerSpacePoints(
      spacePoints, ShowerStartPosition, ShowerDirection);

    //Find the length as the value that contains % of the hits
    int lengthIter = fPercentile * spacePoints.size();

    //Find the length
    double ShowerLength = IShowerTool::GetLArPandoraShowerAlg().SpacePointProjection(
      spacePoints[lengthIter], ShowerStartPosition, ShowerDirection);
    double ShowerMaxProjection = IShowerTool::GetLArPandoraShowerAlg().SpacePointProjection(
      spacePoints[spacePoints.size() - 1], ShowerStartPosition, ShowerDirection);

    double ShowerLengthError = ShowerMaxProjection - ShowerLength;

    //Order the spacepoints in perpendicular
    IShowerTool::GetLArPandoraShowerAlg().OrderShowerSpacePointsPerpendicular(
      spacePoints, ShowerStartPosition, ShowerDirection);

    //Find the length as the value that contains % of the hits
    int perpIter = fPercentile * spacePoints.size();

    //Find the width of the shower
    double ShowerWidth = IShowerTool::GetLArPandoraShowerAlg().SpacePointPerpendicular(
      spacePoints[perpIter], ShowerStartPosition, ShowerDirection);
    // double ShowerMaxWidth = IShowerTool::GetLArPandoraShowerAlg().SpacePointPerpendicular(
    //     spacePoints[spacePoints.size() -1], ShowerStartPosition, ShowerDirection);

    double ShowerAngle = std::atan(ShowerWidth / ShowerLength);
    double ShowerAngleError = -999; //TODO: Do properly

    // Fill the shower element holder
    ShowerEleHolder.SetElement(ShowerLength, ShowerLengthError, fShowerLengthOutputLabel);
    ShowerEleHolder.SetElement(ShowerAngle, ShowerAngleError, fShowerOpeningAngleOutputLabel);

    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerLengthPercentile)
