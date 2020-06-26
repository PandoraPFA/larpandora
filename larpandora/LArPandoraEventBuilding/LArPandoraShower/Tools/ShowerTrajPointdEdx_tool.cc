//############################################################################
//### Name:        ShowerTrajPointdEdx                           ###
//### Author:      Dominic Barker (dominic.barker@sheffield.ac.uk)         ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the dEdx of the start track of the     ###
//###              shower using the standard calomitry module. This        ###
//###              takes the sliding fit trajectory to make a 3D dEdx.     ###
//###              This module is best used with the sliding linear fit    ###
//###              and ShowerTrackTrajToSpacePoint                         ###
//############################################################################

#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft Includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Algs/LArPandoraShowerAlg.h"

//C++ Includes
#include <iostream>
#include <vector>

//Root Includes
#include "TVector3.h"

namespace ShowerRecoTools{

  class ShowerTrajPointdEdx:IShowerTool {

    public:

      ShowerTrajPointdEdx(const fhicl::ParameterSet& pset);

      ~ShowerTrajPointdEdx();

      //Physics Function. Calculate the dEdx.
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder) override;

      void FinddEdxLength(std::vector<double>& dEdx_vec, std::vector<double>& dEdx_val);

    private:

      //Servcies and Algorithms
      art::ServiceHandle<geo::Geometry> fGeom;
      calo::CalorimetryAlg fCalorimetryAlg;

      detinfo::DetectorProperties const* fDetProp;

      //fcl parameters
      float fMinAngleToWire;  //Minimum angle between the wire direction and the shower
      //direction for the spacepoint to be used. Default means
      //the cut has no effect. In radians.
      float fShapingTime;     //Shaping time of the ASIC defualt so we don't cut on track
      //going too much into the plane. In Microseconds
      float fMinDistCutOff;   //Distance in wires a hit has to be from the start position
      //to be used
      float fMaxDist,MaxDist;         //Distance in wires a that a trajectory point can be from a
      //spacepoint to match to it.
      float fdEdxTrackLength,dEdxTrackLength; //Max Distance a spacepoint can be away from the start of the
      //track. In cm
      float fdEdxCut;
      bool fUseMedian;        //Use the median value as the dEdx rather than the mean.
      bool fCutStartPosition; //Remove hits using MinDistCutOff from the vertex as well.

      art::InputTag fPFParticleLabel;
      int           fVerbose;

      std::string fShowerStartPositionInputLabel;
      std::string fInitialTrackSpacePointsInputLabel;
      std::string fInitialTrackInputLabel;
      std::string fShowerdEdxOutputLabel;
      std::string fShowerBestPlaneOutputLabel;
      std::string fShowerdEdxVecOutputLabel;
  };


  ShowerTrajPointdEdx::ShowerTrajPointdEdx(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
    fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>()),
    fMinAngleToWire(pset.get<float>("MinAngleToWire")),
    fShapingTime(pset.get<float>("ShapingTime")),
    fMinDistCutOff(pset.get<float>("MinDistCutOff")),
    fMaxDist(pset.get<float>("MaxDist")),
    fdEdxTrackLength(pset.get<float>("dEdxTrackLength")),
    fdEdxCut(pset.get<float>("dEdxCut")),
    fUseMedian(pset.get<bool>("UseMedian")),
    fCutStartPosition(pset.get<bool>("CutStartPosition")),
    fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel")),
    fVerbose(pset.get<int>("Verbose")),
    fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel")),
    fInitialTrackSpacePointsInputLabel(pset.get<std::string>("InitialTrackSpacePointsInputLabel")),
    fInitialTrackInputLabel(pset.get<std::string>("InitialTrackInputLabel")),
    fShowerdEdxOutputLabel(pset.get<std::string>("ShowerdEdxOutputLabel")),
    fShowerBestPlaneOutputLabel(pset.get<std::string>("ShowerBestPlaneOutputLabel")),
    fShowerdEdxVecOutputLabel(pset.get<std::string>("ShowerdEdxVecOutputLabel"))
  {
  }

  ShowerTrajPointdEdx::~ShowerTrajPointdEdx()
  {
  }

  int ShowerTrajPointdEdx::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event,
      reco::shower::ShowerElementHolder& ShowerEleHolder){

    MaxDist         = fMaxDist;
    dEdxTrackLength = fdEdxTrackLength;

    // Shower dEdx calculation
    if(!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)){
      if (fVerbose)
        mf::LogError("ShowerTrajPointdEdx") << "Start position not set, returning "<< std::endl;
      return 1;
    }
    if(!ShowerEleHolder.CheckElement(fInitialTrackSpacePointsInputLabel)){
      if (fVerbose)
        mf::LogError("ShowerTrajPointdEdx") << "Initial Track Spacepoints is not set returning"<< std::endl;
      return 1;
    }
    if(!ShowerEleHolder.CheckElement(fInitialTrackInputLabel)){
      if (fVerbose)
        mf::LogError("ShowerTrajPointdEdx") << "Initial Track is not set"<< std::endl;
      return 1;
    }

    //Get the initial track hits
    std::vector<art::Ptr<recob::SpacePoint> > tracksps;
    ShowerEleHolder.GetElement(fInitialTrackSpacePointsInputLabel,tracksps);

    if(tracksps.size() == 0){
      if (fVerbose)
        mf::LogWarning("ShowerTrajPointdEdx") << "no spacepointsin the initial track" << std::endl;
      return 0;
    }

    // Get the spacepoints
    art::Handle<std::vector<recob::SpacePoint> > spHandle;
    if (!Event.getByLabel(fPFParticleLabel, spHandle)){
      throw cet::exception("ShowerTrajPointdEdx") << "Could not configure the spacepoint handle. Something is configured incorrectly. Stopping";
      return 1;
    }

    // Get the hits associated with the space points
    art::FindManyP<recob::Hit>& fmsp = ShowerEleHolder.GetFindManyP<recob::Hit>(
        spHandle, Event, fPFParticleLabel);
    if(!fmsp.isValid()){
      throw cet::exception("ShowerTrajPointdEdx") << "Spacepoint and hit association not valid. Stopping.";
      return 1;
    }

    //Only consider hits in the same tpcs as the vertex.
    TVector3 ShowerStartPosition = {-999,-999,-999};
    ShowerEleHolder.GetElement(fShowerStartPositionInputLabel,ShowerStartPosition);
    geo::TPCID vtxTPC = fGeom->FindTPCAtPosition(ShowerStartPosition);

    //Get the initial track
    recob::Track InitialTrack;
    ShowerEleHolder.GetElement(fInitialTrackInputLabel,InitialTrack);

    //Don't care that I could use a vector.
    std::map<int,std::vector<double > > dEdx_vec;
    std::map<int,std::vector<double> >  dEdx_vecErr;
    std::map<int,int> num_hits;

    for(geo::PlaneID plane_id: fGeom->IteratePlaneIDs()){
      dEdx_vec[plane_id.Plane] = {};
      dEdx_vecErr[plane_id.Plane] = {};
      num_hits[plane_id.Plane] = 0;
    }

    //Loop over the spacepoints
    for(auto const sp: tracksps){

      //Get the associated hit
      std::vector<art::Ptr<recob::Hit> > hits = fmsp.at(sp.key());
      if(hits.size() == 0){
        if (fVerbose)
          mf::LogWarning("ShowerTrajPointdEdx") << "no hit for the spacepoint. This suggest the find many is wrong."<< std::endl;
        continue;
      }
      const art::Ptr<recob::Hit> hit = hits[0];
      double wirepitch = fGeom->WirePitch((geo::PlaneID)hit->WireID());

      //Only consider hits in the same tpc
      geo::PlaneID planeid = hit->WireID();
      geo::TPCID TPC = planeid.asTPCID();
      if (TPC !=vtxTPC){continue;}

      //Ignore spacepoints within a few wires of the vertex.
      double dist_from_start = (IShowerTool::GetLArPandoraShowerAlg().SpacePointPosition(sp) - ShowerStartPosition).Mag();

      if(fCutStartPosition){
        if(dist_from_start < fMinDistCutOff*wirepitch){continue;}

        if(dist_from_start > dEdxTrackLength){continue;}
      }

      //Find the closest trajectory point of the track. These should be in order if the user has used ShowerTrackTrajToSpacePoint_tool but the sake of gernicness I'll get the cloest sp.
      unsigned int index = 999;
      double MinDist = 999;
      for(unsigned int traj=0; traj< InitialTrack.NumberTrajectoryPoints(); ++traj){

        geo::Point_t TrajPositionPoint = InitialTrack.LocationAtPoint(traj);
        TVector3 TrajPosition = {TrajPositionPoint.X(),TrajPositionPoint.Y(),TrajPositionPoint.Z()};


        //ignore bogus info.
        auto flags = InitialTrack.FlagsAtPoint(traj);
        if(flags.isSet(recob::TrajectoryPointFlagTraits::NoPoint))
        {continue;}

        TVector3 pos = IShowerTool::GetLArPandoraShowerAlg().SpacePointPosition(sp) - TrajPosition;

        if(pos.Mag() < MinDist && pos.Mag()< MaxDist*wirepitch){
          MinDist = pos.Mag();
          index = traj;
        }
      }

      //If there is no matching trajectory point then bail.
      if(index == 999){continue;}

      geo::Point_t TrajPositionPoint = InitialTrack.LocationAtPoint(index);
      TVector3 TrajPosition = {TrajPositionPoint.X(),TrajPositionPoint.Y(),TrajPositionPoint.Z()};

      geo::Point_t TrajPositionStartPoint = InitialTrack.LocationAtPoint(0);
      TVector3 TrajPositionStart = {TrajPositionStartPoint.X(),TrajPositionStartPoint.Y(),TrajPositionStartPoint.Z()};

      //Ignore values with 0 mag from the start position
      if((TrajPosition - TrajPositionStart).Mag() == 0){continue;}
      if((TrajPosition - ShowerStartPosition).Mag() == 0){continue;}

      if((TrajPosition-TrajPositionStart).Mag() < fMinDistCutOff*wirepitch){continue;}


      //Get the direction of the trajectory point
      geo::Vector_t TrajDirection_vec = InitialTrack.DirectionAtPoint(index);
      TVector3 TrajDirection = {TrajDirection_vec.X(),TrajDirection_vec.Y(),TrajDirection_vec.Z()};

      //If the direction is in the same direction as the wires within some tolerance the hit finding struggles. Let remove these.
      // Note that we project in the YZ plane to make sure we are not cutting on
      // the angle into the wire planes, that should be done by the shaping time cut
      TVector3 TrajDirectionYZ = {0,TrajDirection_vec.Y(),TrajDirection_vec.Z()};
      TVector3 PlaneDirection = fGeom->Plane(planeid).GetIncreasingWireDirection();

      if(TMath::Abs((TMath::Pi()/2 - TrajDirectionYZ.Angle(PlaneDirection))) < fMinAngleToWire){
        if (fVerbose)
          mf::LogWarning("ShowerTrajPointdEdx")
            << "remove from angle cut" << std::endl;
        continue;
      }

      //If the direction is too much into the wire plane then the shaping amplifer cuts the charge. Lets remove these events.
      double velocity = fDetProp->DriftVelocity(fDetProp->Efield(), fDetProp->Temperature());
      double distance_in_x = TrajDirection.X()*(wirepitch/TrajDirection.Dot(PlaneDirection));
      double time_taken = TMath::Abs(distance_in_x/velocity);

      //Shaping time doesn't seem to exist in a global place so add it as a fcl.
      if(fShapingTime < time_taken){
        if (fVerbose)
          mf::LogWarning("ShowerTrajPointdEdx")
            << "move for shaping time" << std::endl;
        continue;
      }

      //Iterate the number of hits on the plane
      ++num_hits[planeid.Plane];

      if((TrajPosition-TrajPositionStart).Mag() > dEdxTrackLength){continue;}


      //If we still exist then we can be used in the calculation. Calculate the 3D pitch
      double trackpitch = (TrajDirection*(wirepitch/TrajDirection.Dot(PlaneDirection))).Mag();

      //Calculate the dQdx
      double dQdx = hit->Integral()/trackpitch;

      //Calculate the dEdx
      double dEdx = fCalorimetryAlg.dEdx_AREA(dQdx, hit->PeakTime(), planeid.Plane);

      //Add the value to the dEdx
      dEdx_vec[planeid.Plane].push_back(dEdx);
    }

    //Search for blow ups and gradient changes.
    //Electrons have a very flat dEdx as function of energy till ~10MeV.
    //If there is a sudden jump particle has probably split
    //If there is very large dEdx we have either calculated it wrong (probably) or the Electron is coming to end.
    //Assumes hits are ordered!
    std::map<int,std::vector<double > > dEdx_vec_cut;

    for(geo::PlaneID plane_id: fGeom->IteratePlaneIDs()){
      dEdx_vec_cut[plane_id.Plane] = {};
    }

    // int max_hits   = -999;
    // int best_plane = -999;
    // for(auto const& dEdx_plane: dEdx_vec){
    //   if((int) dEdx_plane.second.size() > max_hits){
    //     best_plane = dEdx_plane.first;
    //     max_hits   = dEdx_plane.second.size();
    //   }
    // }

    //Choose max hits based on hitnum
    int max_hits   = -999;
    int best_plane = -999;
    for(auto const& num_hits_plane: num_hits){
      if(num_hits_plane.second > max_hits){
        best_plane = num_hits_plane.first;
        max_hits   = num_hits_plane.second;
      }
    }


    for(auto& dEdx_plane: dEdx_vec){
      FinddEdxLength(dEdx_plane.second, dEdx_vec_cut[dEdx_plane.first]);
    }


    //Never have the stats to do a landau fit and get the most probable value. User decides if they want the median value or the mean.
    std::vector<double> dEdx_val;
    std::vector<double> dEdx_valErr;
    for(auto const& dEdx_plane: dEdx_vec_cut){

      if((dEdx_plane.second).size() == 0){
        dEdx_val.push_back(-999);
        dEdx_valErr.push_back(-999);
        continue;
      }

      if(fUseMedian){
        dEdx_val.push_back(TMath::Median((dEdx_plane.second).size(), &(dEdx_plane.second)[0]));
      }
      else{
        //Else calculate the mean value.
        double dEdx_mean = 0;
        for(auto const& dEdx: dEdx_plane.second){
          if(dEdx > 10 || dEdx < 0){continue;}
          dEdx_mean += dEdx;
        }
        dEdx_val.push_back(dEdx_mean/(float)(dEdx_plane.second).size());
      }
    }

    if (fVerbose>1){
      for(auto const& plane: dEdx_vec_cut){
        std::cout << "#Plane: " << plane.first << " #" << std::endl;
        for(auto const& dEdx: plane.second){
          std::cout << "dEdx: " << dEdx << std::endl;
        }
      }
    }

    //Need to sort out errors sensibly.
    ShowerEleHolder.SetElement(dEdx_val,dEdx_valErr,fShowerdEdxOutputLabel);
    ShowerEleHolder.SetElement(best_plane,fShowerBestPlaneOutputLabel);
    ShowerEleHolder.SetElement(dEdx_vec_cut,fShowerdEdxVecOutputLabel);
    return 0;
  }


  void ShowerTrajPointdEdx::FinddEdxLength(std::vector<double>& dEdx_vec, std::vector<double>& dEdx_val){

    //As default do not apply this cut.
    if(fdEdxCut > 10){
      dEdx_val = dEdx_vec;
      return;
    }

    //Can only do this with 4 hits.
    if(dEdx_vec.size() < 4){
      dEdx_val = dEdx_vec;
      return;
    }

    bool upperbound = false;

    //See if we are in the upper bound or upper bound defined by the cut.
    int upperbound_int = 0;
    if(dEdx_vec[0] > fdEdxCut){++upperbound_int;}
    if(dEdx_vec[1] > fdEdxCut){++upperbound_int;}
    if(dEdx_vec[2] > fdEdxCut){++upperbound_int;}
    if(upperbound_int > 1){upperbound = true;}


    dEdx_val.push_back(dEdx_vec[0]);
    dEdx_val.push_back(dEdx_vec[1]);
    dEdx_val.push_back(dEdx_vec[2]);

    for(unsigned int dEdx_iter=2; dEdx_iter<dEdx_vec.size(); ++dEdx_iter){

      //The Function of dEdx as a function of E is flat above ~10 MeV.
      //We are looking for a jump up (or down) above the ladau width in the dEx
      //to account account for pair production.
      //Dom Estimates that the somwhere above 0.28 MeV will be a good cut but 999 will prevent this stage.
      double dEdx = dEdx_vec[dEdx_iter];

      //We are really poo at physics and so attempt to find the pair production
      if(upperbound){
        if(dEdx > fdEdxCut){
          dEdx_val.push_back(dEdx);
          if (fVerbose>1)
            std::cout << "Adding dEdx: "<< dEdx<< std::endl;
          continue;
        }
        else{
          //Maybe its a landau fluctation lets try again.
          if(dEdx_iter < dEdx_vec.size()-1){
            if(dEdx_vec[dEdx_iter+1] > fdEdxCut){
              if (fVerbose>1)
                std::cout << "Next dEdx hit is good removing hit"<< dEdx<< std::endl;
              continue;
            }
          }
          //I'll let one more value
          if(dEdx_iter<dEdx_vec.size()-2){
            if(dEdx_vec[dEdx_iter+2] > fdEdxCut){
              if (fVerbose>1)
                std::cout << "Next Next dEdx hit is good removing hit"<< dEdx<< std::endl;
              continue;
            }
          }
          //We are hopefully we have one of our electrons has died.
          break;
        }
      }
      else{
        if(dEdx < fdEdxCut){
          dEdx_val.push_back(dEdx);
          if (fVerbose>1)
            std::cout << "Adding dEdx: "<< dEdx<< std::endl;
          continue;
        }
        else{
          //Maybe its a landau fluctation lets try again.
          if(dEdx_iter < dEdx_vec.size()-1){
            if(dEdx_vec[dEdx_iter+1] > fdEdxCut){
              if (fVerbose>1)
                std::cout << "Next dEdx hit is good removing hit "<< dEdx<< std::endl;
              continue;
            }
          }
          //I'll let one more value
          if(dEdx_iter < dEdx_vec.size()-2){
            if(dEdx_vec[dEdx_iter+2] > fdEdxCut){
              if (fVerbose>1)
                std::cout << "Next Next dEdx hit is good removing hit "<< dEdx<< std::endl;
              continue;
            }
          }
          //We are hopefully in the the pair production zone.
          break;
        }
      }
    }
    return;
  }

}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerTrajPointdEdx)

