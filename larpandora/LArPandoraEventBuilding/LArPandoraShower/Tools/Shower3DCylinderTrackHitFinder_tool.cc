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
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Algs/LArPandoraShowerAlg.h"

namespace ShowerRecoTools{

  class Shower3DCylinderTrackHitFinder:IShowerTool {
    public:

      Shower3DCylinderTrackHitFinder(const fhicl::ParameterSet& pset);

      ~Shower3DCylinderTrackHitFinder();

      //Generic Track Finder
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder
          ) override;

    private:

      std::vector<art::Ptr<recob::SpacePoint> > FindTrackSpacePoints(std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
          TVector3& showerStartPosition,TVector3& showerDirection);


      //Fcl paramters
      float fMaxProjectionDist;    //Maximum projection along shower direction.
      float fMaxPerpendicularDist; //Maximum perpendicular distance, radius of cylinder
      bool  fForwardHitsOnly;      //Only take hits downstream of shower vertex
      //(projection>0)

      art::InputTag fPFParticleLabel;
      int           fVerbose;

      std::string fShowerStartPositionInputLabel;
      std::string fInitialTrackHitsOutputLabel;
      std::string fInitialTrackSpacePointsOutputLabel;
      std::string fShowerDirectionInputLabel;
  };


  Shower3DCylinderTrackHitFinder::Shower3DCylinderTrackHitFinder(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fMaxProjectionDist(pset.get<float>("MaxProjectionDist")),
    fMaxPerpendicularDist(pset.get<float>("MaxPerpendicularDist")),
    fForwardHitsOnly(pset.get<bool>("ForwardHitsOnly")),
    fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel")),
    fVerbose(pset.get<int>("Verbose")),
    fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel")),
    fInitialTrackHitsOutputLabel(pset.get<std::string>("InitialTrackHitsOutputLabel")),
    fInitialTrackSpacePointsOutputLabel(pset.get<std::string>("InitialTrackSpacePointsOutputLabel")),
    fShowerDirectionInputLabel(pset.get<std::string>("ShowerDirectionInputLabel"))
  {
  }

  Shower3DCylinderTrackHitFinder::~Shower3DCylinderTrackHitFinder()
  {
  }

  int Shower3DCylinderTrackHitFinder::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder){

    //This is all based on the shower vertex being known. If it is not lets not do the track
    if(!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)){
      if (fVerbose)
        mf::LogError("Shower3DCylinderTrackHitFinder") << "Start position not set, returning "<< std::endl;
      return 1;
    }
    if(!ShowerEleHolder.CheckElement("ShowerDirection")){
      if (fVerbose)
        mf::LogError("Shower3DCylinderTrackHitFinder") << "Direction not set, returning "<< std::endl;
      return 1;
    }

    TVector3 ShowerStartPosition = {-999,-999,-999};
    ShowerEleHolder.GetElement(fShowerStartPositionInputLabel,ShowerStartPosition);

    TVector3 ShowerDirection     = {-999,-999,-999};
    ShowerEleHolder.GetElement(fShowerDirectionInputLabel,ShowerDirection);

    // Get the assocated pfParicle Handle
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleLabel, pfpHandle)){
      throw cet::exception("Shower3DCylinderTrackHitFinder") << "Could not get the pandora pf particles. Something is not cofingured correctly Please give the correct pandoa module label. Stopping";
      return 1;
    }

    // Get the spacepoint - PFParticle assn
    art::FindManyP<recob::SpacePoint>& fmspp = ShowerEleHolder.GetFindManyP<recob::SpacePoint>(
        pfpHandle, Event, fPFParticleLabel);
    if (!fmspp.isValid()){
      throw cet::exception("Shower3DCylinderTrackHitFinder") << "Trying to get the spacepoint and failed. Something is not configured correctly. Stopping ";
      return 1;
    }

    // Get the spacepoints
    art::Handle<std::vector<recob::SpacePoint> > spHandle;
    if (!Event.getByLabel(fPFParticleLabel, spHandle)){
      throw cet::exception("Shower3DCylinderTrackHitFinder") << "Could not configure the spacepoint handle. Something is configured incorrectly. Stopping";
      return 1;
    }

    // Get the hits associated with the space points
    art::FindManyP<recob::Hit>& fmhsp = ShowerEleHolder.GetFindManyP<recob::Hit>(
        spHandle, Event, fPFParticleLabel);

    // art::FindOneP<recob::Hit> fohsp(spHandle, Event, fPFParticleLabel);
    if(!fmhsp.isValid()){
      throw cet::exception("Shower3DCylinderTrackHitFinder") << "Spacepoint and hit association not valid. Stopping.";
      return 1;
    }

    // Get the SpacePoints
    std::vector<art::Ptr<recob::SpacePoint> > spacePoints = fmspp.at(pfparticle.key());

    //We cannot progress with no spacepoints.
    if(spacePoints.size() == 0){
      if (fVerbose)
        mf::LogError("Shower3DCylinderTrackHitFinder") << "No space points, returning "<< std::endl;
      return 1;
    }

    // Order the spacepoints
    IShowerTool::GetLArPandoraShowerAlg().OrderShowerSpacePoints(spacePoints,ShowerStartPosition,ShowerDirection);

    // Get only the space points from the track
    std::vector<art::Ptr<recob::SpacePoint> > trackSpacePoints;
    trackSpacePoints = FindTrackSpacePoints(spacePoints,ShowerStartPosition,ShowerDirection);

    // Get the hits associated to the space points and seperate them by planes
    std::vector<art::Ptr<recob::Hit> > trackHits;
    for(auto const& spacePoint: trackSpacePoints){
      const art::Ptr<recob::Hit> hit = fmhsp.at(spacePoint.key()).front();
      // const art::Ptr<recob::Hit> hit = fohsp.at(spacePoint.key());
      trackHits.push_back(hit);
    }

    ShowerEleHolder.SetElement(trackHits, fInitialTrackHitsOutputLabel);
    ShowerEleHolder.SetElement(trackSpacePoints,fInitialTrackSpacePointsOutputLabel);

    return 0;
  }

  std::vector<art::Ptr<recob::SpacePoint> > Shower3DCylinderTrackHitFinder::FindTrackSpacePoints(
      std::vector<art::Ptr<recob::SpacePoint> >& spacePoints, TVector3& showerStartPosition,
      TVector3& showerDirection){

    // Make a vector to hold the output space points
    std::vector<art::Ptr<recob::SpacePoint> > trackSpacePoints;

    for (const auto& spacePoint : spacePoints){
      // Calculate the projection along direction and perpendicular distance
      // from "axis" of shower TODO: change alg to return a pair for efficiency
      double proj = IShowerTool::GetLArPandoraShowerAlg().SpacePointProjection(spacePoint,
          showerStartPosition, showerDirection);
      double perp = IShowerTool::GetLArPandoraShowerAlg().SpacePointPerpendicular(spacePoint,
          showerStartPosition, showerDirection, proj);

      if (fForwardHitsOnly && proj<0)
        continue;

      if (TMath::Abs(proj)<fMaxProjectionDist && TMath::Abs(perp)<fMaxPerpendicularDist)
        trackSpacePoints.push_back(spacePoint);
    }
    return trackSpacePoints;
  }

}


DEFINE_ART_CLASS_TOOL(ShowerRecoTools::Shower3DCylinderTrackHitFinder)