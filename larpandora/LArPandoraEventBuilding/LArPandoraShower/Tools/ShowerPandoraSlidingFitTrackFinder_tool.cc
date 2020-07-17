//#############################################################################
//### Name:        ShowerPandoraSlidingFitTrackFinder                       ###
//### Author:      Dominic Barker (dominic.barker@sheffield.ac.uk           ###
//### Date:        30.07.19                                                 ###
//### Description: Tool for finding the initial shower track using the      ###
//###              pandora sliding fit calculation. This method is derived  ###
//###              from the PandoraTrackCreationModule.cc                   ###
//#############################################################################

//Framework Includes
#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"

namespace ShowerRecoTools{

  class ShowerPandoraSlidingFitTrackFinder:IShowerTool {
    public:

      ShowerPandoraSlidingFitTrackFinder(const fhicl::ParameterSet& pset);

      ~ShowerPandoraSlidingFitTrackFinder();

      //Generic Track Finder
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder) override;
    private:

      void InitialiseProducers() override;

      //Function to add the assoctions
      int AddAssociations(const art::Ptr<recob::PFParticle>& pfpPtr, art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder) override;


      // Define standard art tool interface.
      art::ServiceHandle<geo::Geometry> fGeom;

      //fcl paramaters
      int   fVerbose;
      float fSlidingFitHalfWindow; //To Describe
      float fMinTrajectoryPoints;  //Minimum number of trajectory point to say the track is good.
      std::string fInitialTrackOutputLabel;
      std::string fInitialTrackLengthOutputLabel;
      std::string fShowerStartPositionInputLabel;
      std::string fShowerDirectionInputLabel;
      std::string fInitialTrackSpacePointsInputLabel;
      std::string fInitialTrackHitsInputLabel;
  };


  ShowerPandoraSlidingFitTrackFinder::ShowerPandoraSlidingFitTrackFinder(const fhicl::ParameterSet& pset):
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fVerbose(pset.get<int>("Verbose")),
    fSlidingFitHalfWindow(pset.get<float> ("SlidingFitHalfWindow")),
    fMinTrajectoryPoints(pset.get<float> ("MinTrajectoryPoints")),
    fInitialTrackOutputLabel(pset.get<std::string>("InitialTrackOutputLabel")),
    fInitialTrackLengthOutputLabel(pset.get<std::string>("InitialTrackLengthOutputLabel")),
    fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel")),
    fShowerDirectionInputLabel(pset.get<std::string>("ShowerDirectionInputLabel")),
    fInitialTrackSpacePointsInputLabel(pset.get<std::string>("InitialTrackSpacePointsInputLabel")),
    fInitialTrackHitsInputLabel(pset.get<std::string>("InitialTrackHitsInputLabel"))
  {
  }

  ShowerPandoraSlidingFitTrackFinder::~ShowerPandoraSlidingFitTrackFinder()
  {
  }

  void ShowerPandoraSlidingFitTrackFinder::InitialiseProducers(){

    InitialiseProduct<std::vector<recob::Track> >(fInitialTrackOutputLabel);
    InitialiseProduct<art::Assns<recob::Shower, recob::Track > >("ShowerTrackAssn");
    InitialiseProduct<art::Assns<recob::Track, recob::Hit > >("ShowerTrackHitAssn");

  }


  //This whole idea is stolen from PandoraTrackCreationModule so credit goes to the Pandora guys.
  int ShowerPandoraSlidingFitTrackFinder::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event,
      reco::shower::ShowerElementHolder& ShowerEleHolder
      ){
    //This is all based on the shower vertex being known. If it is not lets not do the track
    if(!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)){
      if (fVerbose)
        mf::LogError("ShowerPandoraSlidingFitTrackFinder") << "Start position not set, returning "<< std::endl;
      return 1;
    }
    if(!ShowerEleHolder.CheckElement(fShowerDirectionInputLabel)){
      if (fVerbose)
        mf::LogError("ShowerPandoraSlidingFitTrackFinder") << "Direction not set, returning "<< std::endl;
      return 1;
    }
    if(!ShowerEleHolder.CheckElement(fInitialTrackSpacePointsInputLabel)){
      if (fVerbose)
        mf::LogError("ShowerPandoraSlidingFitTrackFinder") << "Initial Spacepoints not set, returning "<< std::endl;
      return 1;
    }

    TVector3 ShowerStartPosition = {-999,-999,-999};
    ShowerEleHolder.GetElement(fShowerStartPositionInputLabel,ShowerStartPosition);

    TVector3 ShowerDirection     = {-999,-999,-999};
    ShowerEleHolder.GetElement(fShowerDirectionInputLabel,ShowerDirection);

    std::vector<art::Ptr<recob::SpacePoint> > spacepoints;
    ShowerEleHolder.GetElement(fInitialTrackSpacePointsInputLabel,spacepoints);

    const unsigned int nWirePlanes(fGeom->MaxPlanes());

    if (nWirePlanes > 3)
      throw cet::exception("LArPandoraTrackCreation") << " LArPandoraTrackCreation::produce --- More than three wire planes present ";

    if ((0 == fGeom->Ncryostats()) || (0 == fGeom->NTPC(0)))
      throw cet::exception("LArPandoraTrackCreation") << " LArPandoraTrackCreation::produce --- unable to access first tpc in first cryostat ";

    std::unordered_set<geo::_plane_proj> planeSet;
    for (unsigned int iPlane = 0; iPlane < nWirePlanes; ++iPlane)
      (void) planeSet.insert(fGeom->TPC(0, 0).Plane(iPlane).View());

    // ATTN: Expectations here are that the input geometry corresponds to either a single or dual phase LArTPC.  For single phase we expect
    // three views, U, V and either W or Y, for dual phase we expect two views, W and Y.
    const bool isDualPhase(fGeom->MaxPlanes() == 2);

    if (nWirePlanes != planeSet.size())
      throw cet::exception("LArPandoraTrackCreation") << " LArPandoraGeometry::LoadGeometry --- geometry description for wire plane(s) missing ";

    if (isDualPhase && (!planeSet.count(geo::kW) || !planeSet.count(geo::kY)))
      throw cet::exception("LArPandoraTrackCreation") << " LArPandoraGeometry::LoadGeometry --- dual phase scenario; expect to find w and y views ";

    if (!isDualPhase && (!planeSet.count(geo::kU) || !planeSet.count(geo::kV) || (planeSet.count(geo::kW) && planeSet.count(geo::kY))))
      throw cet::exception("LArPandoraTrackCreation") << " LArPandoraGeometry::LoadGeometry --- single phase scenatio; expect to find u and v views; if there is one further view, it must be w or y ";

    const bool useYPlane((nWirePlanes > 2) && planeSet.count(geo::kY));

    // ATTN: In the dual phase mode, map the wire planes as follows W->U and Y->V.  This mapping was chosen so that the dual phase wire
    // planes, which are inherently induction only, are mapped to induction planes in the single phase geometry.
    const float wirePitchU(fGeom->WirePitch((isDualPhase ? geo::kW : geo::kU)));
    const float wirePitchV(fGeom->WirePitch((isDualPhase ? geo::kY : geo::kV)));
    const float wirePitchW((nWirePlanes < 3) ? 0.5f * (wirePitchU + wirePitchV) : (useYPlane) ? fGeom->WirePitch(geo::kY) :
        fGeom->WirePitch(geo::kW));


    const pandora::CartesianVector vertexPosition(ShowerStartPosition.X(), ShowerStartPosition.Y(),
        ShowerStartPosition.Z());

    pandora::CartesianPointVector cartesianPointVector;
    for (const art::Ptr<recob::SpacePoint> spacePoint : spacepoints)
      cartesianPointVector.emplace_back(pandora::CartesianVector(spacePoint->XYZ()[0],
            spacePoint->XYZ()[1], spacePoint->XYZ()[2]));

    lar_content::LArTrackStateVector trackStateVector;
    pandora::IntVector indexVector;
    try{
      lar_content::LArPfoHelper::GetSlidingFitTrajectory(cartesianPointVector, vertexPosition,
          fSlidingFitHalfWindow, wirePitchW, trackStateVector, &indexVector);
    }
    catch (const pandora::StatusCodeException &){
      if (fVerbose)
        mf::LogWarning("ShowerPandoraSlidingFitTrackFinder") << "Unable to extract sliding fit trajectory" << std::endl;
      return 1;
    }
    if (trackStateVector.size() < fMinTrajectoryPoints){
      if (fVerbose)
        mf::LogWarning("ShowerPandoraSlidingFitTrackFinder") << "Insufficient input trajectory points to build track: " << trackStateVector.size();
      return 1;
    }

    if (trackStateVector.empty())
      throw cet::exception("ShowerPandoraSlidingFitTrackFinder") << "BuildTrack - No input trajectory points provided " << std::endl;

    recob::tracking::Positions_t xyz;
    recob::tracking::Momenta_t pxpypz;
    recob::TrackTrajectory::Flags_t flags;

    for (const lar_content::LArTrackState &trackState : trackStateVector)
    {
      xyz.emplace_back(recob::tracking::Point_t(trackState.GetPosition().GetX(),
            trackState.GetPosition().GetY(), trackState.GetPosition().GetZ()));
      pxpypz.emplace_back(recob::tracking::Vector_t(trackState.GetDirection().GetX(),
            trackState.GetDirection().GetY(), trackState.GetDirection().GetZ()));

      // Set flag NoPoint if point has bogus coordinates, otherwise use clean flag set
      if (std::fabs(trackState.GetPosition().GetX()-util::kBogusF)<std::numeric_limits<float>::epsilon() &&
          std::fabs(trackState.GetPosition().GetY()-util::kBogusF)<std::numeric_limits<float>::epsilon() &&
          std::fabs(trackState.GetPosition().GetZ()-util::kBogusF)<std::numeric_limits<float>::epsilon())
      {
        flags.emplace_back(recob::TrajectoryPointFlags(recob::TrajectoryPointFlags::InvalidHitIndex,
              recob::TrajectoryPointFlagTraits::NoPoint));
      } else {
        flags.emplace_back(recob::TrajectoryPointFlags());
      }
    }

    // note from gc: eventually we should produce a TrackTrajectory, not a Track with empty covariance matrix and bogus chi2, etc.
    recob::Track InitialTrack  =  recob::Track(recob::TrackTrajectory(std::move(xyz),std::move(pxpypz), std::move(flags), false),
        util::kBogusI, util::kBogusF, util::kBogusI, recob::tracking::SMatrixSym55(),
        recob::tracking::SMatrixSym55(), pfparticle.key());

    ShowerEleHolder.SetElement(InitialTrack,fInitialTrackOutputLabel);

    TVector3 Start = {InitialTrack.Start().X(), InitialTrack.Start().Y(), InitialTrack.Start().Z()};
    TVector3 End   = {InitialTrack.End().X(), InitialTrack.End().Y(),InitialTrack.End().Z()};
    float tracklength = (Start-End).Mag();

    ShowerEleHolder.SetElement(tracklength,fInitialTrackLengthOutputLabel);

    return 0;
  }


  int ShowerPandoraSlidingFitTrackFinder::AddAssociations(const art::Ptr<recob::PFParticle>& pfpPtr, art::Event& Event,
      reco::shower::ShowerElementHolder& ShowerEleHolder
      ){

    //Check the track has been set
    if(!ShowerEleHolder.CheckElement(fInitialTrackOutputLabel)){
      if (fVerbose)
        mf::LogError("ShowerPandoraSlidingFitTrackFinderAddAssn") << "Track not set so the assocation can not be made  "<< std::endl;
      return 1;
    }

    //Get the size of the ptr as it is.
    int trackptrsize = GetVectorPtrSize(fInitialTrackOutputLabel);

    const art::Ptr<recob::Track> trackptr = GetProducedElementPtr<recob::Track>(fInitialTrackOutputLabel,
        ShowerEleHolder,trackptrsize-1);
    const art::Ptr<recob::Shower> showerptr = GetProducedElementPtr<recob::Shower>("shower",
        ShowerEleHolder);

    AddSingle<art::Assns<recob::Shower, recob::Track> >(showerptr,trackptr,"ShowerTrackAssn");

    std::vector<art::Ptr<recob::Hit> > TrackHits;
    ShowerEleHolder.GetElement(fInitialTrackHitsInputLabel,TrackHits);

    for(auto const& TrackHit: TrackHits){
      AddSingle<art::Assns<recob::Track, recob::Hit> >(trackptr,TrackHit,"ShowerTrackHitAssn");
    }

    return 0;
  }

}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerPandoraSlidingFitTrackFinder)
