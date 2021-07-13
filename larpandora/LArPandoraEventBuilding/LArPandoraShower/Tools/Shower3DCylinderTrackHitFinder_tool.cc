//############################################################################
//### Name:        Shower3DCylinderTrackHitFinder                          ###
//### Author:      Ed Tyley                                                ###
//### Date:        14.06.19                                                ###
//### Description: Tool for finding the initial shower track using 3D      ###
//###              spacepoints within a cylinder along the shower          ###
//###              direction. fcl parameters define cylinder dimensions    ###
//############################################################################

//Framework Includes
#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Algs/LArPandoraShowerAlg.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"

namespace ShowerRecoTools {

  class Shower3DCylinderTrackHitFinder : IShowerTool {
  public:
    Shower3DCylinderTrackHitFinder(const fhicl::ParameterSet& pset);

    //Generic Track Finder
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerEleHolder) override;

  private:
    std::vector<art::Ptr<recob::SpacePoint>> FindTrackSpacePoints(
      std::vector<art::Ptr<recob::SpacePoint>>& spacePoints,
      TVector3& showerStartPosition,
      TVector3& showerDirection);

    //Fcl paramters
    float fMaxProjectionDist;    //Maximum projection along shower direction.
    float fMaxPerpendicularDist; //Maximum perpendicular distance, radius of cylinder
    bool fForwardHitsOnly;       //Only take hits downstream of shower vertex
    //(projection>0)

    art::InputTag fPFParticleLabel;
    int fVerbose;

    std::string fShowerStartPositionInputLabel;
    std::string fInitialTrackHitsOutputLabel;
    std::string fInitialTrackSpacePointsOutputLabel;
    std::string fShowerDirectionInputLabel;
  };

  Shower3DCylinderTrackHitFinder::Shower3DCylinderTrackHitFinder(const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fMaxProjectionDist(pset.get<float>("MaxProjectionDist"))
    , fMaxPerpendicularDist(pset.get<float>("MaxPerpendicularDist"))
    , fForwardHitsOnly(pset.get<bool>("ForwardHitsOnly"))
    , fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel"))
    , fVerbose(pset.get<int>("Verbose"))
    , fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel"))
    , fInitialTrackHitsOutputLabel(pset.get<std::string>("InitialTrackHitsOutputLabel"))
    , fInitialTrackSpacePointsOutputLabel(
        pset.get<std::string>("InitialTrackSpacePointsOutputLabel"))
    , fShowerDirectionInputLabel(pset.get<std::string>("ShowerDirectionInputLabel"))
  {}

  int
  Shower3DCylinderTrackHitFinder::CalculateElement(
    const art::Ptr<recob::PFParticle>& pfparticle,
    art::Event& Event,
    reco::shower::ShowerElementHolder& ShowerEleHolder)
  {

    //This is all based on the shower vertex being known. If it is not lets not do the track
    if (!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)) {
      if (fVerbose)
        mf::LogError("Shower3DCylinderTrackHitFinder")
          << "Start position not set, returning " << std::endl;
      return 1;
    }
    if (!ShowerEleHolder.CheckElement("ShowerDirection")) {
      if (fVerbose)
        mf::LogError("Shower3DCylinderTrackHitFinder")
          << "Direction not set, returning " << std::endl;
      return 1;
    }

    TVector3 ShowerStartPosition = {-999, -999, -999};
    ShowerEleHolder.GetElement(fShowerStartPositionInputLabel, ShowerStartPosition);

    TVector3 ShowerDirection = {-999, -999, -999};
    ShowerEleHolder.GetElement(fShowerDirectionInputLabel, ShowerDirection);

    // Get the assocated pfParicle Handle
    auto const pfpHandle = Event.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);

    // Get the spacepoint - PFParticle assn
    const art::FindManyP<recob::SpacePoint>& fmspp =
      ShowerEleHolder.GetFindManyP<recob::SpacePoint>(pfpHandle, Event, fPFParticleLabel);

    // Get the spacepoints
    auto const spHandle = Event.getValidHandle<std::vector<recob::SpacePoint>>(fPFParticleLabel);

    // Get the hits associated with the space points
    const art::FindManyP<recob::Hit>& fmhsp =
      ShowerEleHolder.GetFindManyP<recob::Hit>(spHandle, Event, fPFParticleLabel);

    // Get the SpacePoints
    std::vector<art::Ptr<recob::SpacePoint>> spacePoints = fmspp.at(pfparticle.key());

    //We cannot progress with no spacepoints.
    if (spacePoints.empty()) {
      if (fVerbose)
        mf::LogError("Shower3DCylinderTrackHitFinder")
          << "No space points, returning " << std::endl;
      return 1;
    }

    // Order the spacepoints
    IShowerTool::GetLArPandoraShowerAlg().OrderShowerSpacePoints(
      spacePoints, ShowerStartPosition, ShowerDirection);

    // Get only the space points from the track
    std::vector<art::Ptr<recob::SpacePoint>> trackSpacePoints;
    trackSpacePoints = FindTrackSpacePoints(spacePoints, ShowerStartPosition, ShowerDirection);

    // Get the hits associated to the space points and seperate them by planes
    std::vector<art::Ptr<recob::Hit>> trackHits;
    for (auto const& spacePoint : trackSpacePoints) {
      const art::Ptr<recob::Hit> hit = fmhsp.at(spacePoint.key()).front();
      // const art::Ptr<recob::Hit> hit = fohsp.at(spacePoint.key());
      trackHits.push_back(hit);
    }

    ShowerEleHolder.SetElement(trackHits, fInitialTrackHitsOutputLabel);
    ShowerEleHolder.SetElement(trackSpacePoints, fInitialTrackSpacePointsOutputLabel);

    return 0;
  }

  std::vector<art::Ptr<recob::SpacePoint>>
  Shower3DCylinderTrackHitFinder::FindTrackSpacePoints(
    std::vector<art::Ptr<recob::SpacePoint>>& spacePoints,
    TVector3& showerStartPosition,
    TVector3& showerDirection)
  {

    // Make a vector to hold the output space points
    std::vector<art::Ptr<recob::SpacePoint>> trackSpacePoints;

    for (const auto& spacePoint : spacePoints) {
      // Calculate the projection along direction and perpendicular distance
      // from "axis" of shower TODO: change alg to return a pair for efficiency
      double proj = IShowerTool::GetLArPandoraShowerAlg().SpacePointProjection(
        spacePoint, showerStartPosition, showerDirection);
      double perp = IShowerTool::GetLArPandoraShowerAlg().SpacePointPerpendicular(
        spacePoint, showerStartPosition, showerDirection, proj);

      if (fForwardHitsOnly && proj < 0) continue;

      if (std::abs(proj) < fMaxProjectionDist && std::abs(perp) < fMaxPerpendicularDist)
        trackSpacePoints.push_back(spacePoint);
    }
    return trackSpacePoints;
  }

}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::Shower3DCylinderTrackHitFinder)
