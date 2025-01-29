//############################################################################
//### Name:        ShowerTrajPointdEdx                                     ###
//### Author:      Dominic Barker (dominic.barker@sheffield.ac.uk)         ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the dEdx of the start track of the     ###
//###              shower using the standard calomitry module. This        ###
//###              takes the sliding fit trajectory to make a 3D dEdx.     ###
//###              This module is best used with the sliding linear fit    ###
//###              and ShowerTrackTrajToSpacePoint                         ###
//############################################################################

//Framework Includes
#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

//ROOT
#include "Math/VectorUtil.h"

using ROOT::Math::VectorUtil::Angle;

namespace ShowerRecoTools {

  class ShowerTrajPointdEdx : IShowerTool {

  public:
    ShowerTrajPointdEdx(const fhicl::ParameterSet& pset);

    //Physics Function. Calculate the dEdx.
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerEleHolder) override;

    void FinddEdxLength(std::vector<double>& dEdx_vec, std::vector<double>& dEdx_val);

  private:
    //Servcies and Algorithms
    art::ServiceHandle<geo::Geometry> fGeom;
    geo::WireReadoutGeom const& fChannelMap = art::ServiceHandle<geo::WireReadout>()->Get();
    calo::CalorimetryAlg fCalorimetryAlg;

    //fcl parameters
    float fMinAngleToWire; //Minimum angle between the wire direction and the shower
    //direction for the spacepoint to be used. Default means
    //the cut has no effect. In radians.
    float fShapingTime; //Shaping time of the ASIC defualt so we don't cut on track
    //going too much into the plane. In Microseconds
    float fMinDistCutOff; //Distance in wires a hit has to be from the start position
    //to be used
    float fMaxDist, MaxDist; //Distance in wires a that a trajectory point can be from a
    //spacepoint to match to it.
    float fdEdxTrackLength,
      dEdxTrackLength; //Max Distance a spacepoint can be away from the start of the
    //track. In cm
    float fdEdxCut;
    bool fUseMedian;        //Use the median value as the dEdx rather than the mean.
    bool fCutStartPosition; //Remove hits using MinDistCutOff from the vertex as well.

    bool fT0Correct;       // Whether to look for a T0 associated to the PFP
    bool fSCECorrectPitch; // Whether to correct the "squeezing" of pitch, requires corrected input
    bool
      fSCECorrectEField; // Whether to use the local electric field, from SpaceChargeService, in recombination calc.
    bool
      fSCEInputCorrected; // Whether the input has already been corrected for spatial SCE distortions

    bool fSumHitSnippets; // Whether to treat hits individually or only one hit per snippet

    art::InputTag fPFParticleLabel;
    int fVerbose;

    std::string fShowerStartPositionInputLabel;
    std::string fInitialTrackHitsInputLabel;
    std::string fInitialTrackSpacePointsInputLabel;
    std::string fInitialTrackInputLabel;
    std::string fShowerdEdxOutputLabel;
    std::string fShowerBestPlaneOutputLabel;
    std::string fShowerdEdxVecOutputLabel;
  };

  ShowerTrajPointdEdx::ShowerTrajPointdEdx(const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
    , fMinAngleToWire(pset.get<float>("MinAngleToWire"))
    , fShapingTime(pset.get<float>("ShapingTime"))
    , fMinDistCutOff(pset.get<float>("MinDistCutOff"))
    , fMaxDist(pset.get<float>("MaxDist"))
    , fdEdxTrackLength(pset.get<float>("dEdxTrackLength"))
    , fdEdxCut(pset.get<float>("dEdxCut"))
    , fUseMedian(pset.get<bool>("UseMedian"))
    , fCutStartPosition(pset.get<bool>("CutStartPosition"))
    , fT0Correct(pset.get<bool>("T0Correct"))
    , fSCECorrectPitch(pset.get<bool>("SCECorrectPitch"))
    , fSCECorrectEField(pset.get<bool>("SCECorrectEField"))
    , fSCEInputCorrected(pset.get<bool>("SCEInputCorrected"))
    , fSumHitSnippets(pset.get<bool>("SumHitSnippets"))
    , fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel"))
    , fVerbose(pset.get<int>("Verbose"))
    , fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel"))
    , fInitialTrackHitsInputLabel(pset.get<std::string>("InitialTrackHitsInputLabel"))
    , fInitialTrackSpacePointsInputLabel(pset.get<std::string>("InitialTrackSpacePointsInputLabel"))
    , fInitialTrackInputLabel(pset.get<std::string>("InitialTrackInputLabel"))
    , fShowerdEdxOutputLabel(pset.get<std::string>("ShowerdEdxOutputLabel"))
    , fShowerBestPlaneOutputLabel(pset.get<std::string>("ShowerBestPlaneOutputLabel"))
    , fShowerdEdxVecOutputLabel(pset.get<std::string>("ShowerdEdxVecOutputLabel"))
  {
    if ((fSCECorrectPitch || fSCECorrectEField) && !fSCEInputCorrected) {
      throw cet::exception("ShowerTrajPointdEdx")
        << "Can only correct for SCE if input is already corrected" << std::endl;
    }
  }

  int ShowerTrajPointdEdx::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                                            art::Event& Event,
                                            reco::shower::ShowerElementHolder& ShowerEleHolder)
  {

    MaxDist = fMaxDist;
    dEdxTrackLength = fdEdxTrackLength;

    // Shower dEdx calculation
    if (!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)) {
      if (fVerbose)
        mf::LogError("ShowerTrajPointdEdx") << "Start position not set, returning " << std::endl;
      return 1;
    }
    if (!ShowerEleHolder.CheckElement(fInitialTrackSpacePointsInputLabel)) {
      if (fVerbose)
        mf::LogError("ShowerTrajPointdEdx")
          << "Initial Track Spacepoints is not set returning" << std::endl;
      return 1;
    }
    if (!ShowerEleHolder.CheckElement(fInitialTrackInputLabel)) {
      if (fVerbose) mf::LogError("ShowerTrajPointdEdx") << "Initial Track is not set" << std::endl;
      return 1;
    }

    //Get the initial track hits
    std::vector<art::Ptr<recob::SpacePoint>> tracksps;
    ShowerEleHolder.GetElement(fInitialTrackSpacePointsInputLabel, tracksps);

    if (tracksps.empty()) {
      if (fVerbose)
        mf::LogWarning("ShowerTrajPointdEdx") << "no spacepointsin the initial track" << std::endl;
      return 0;
    }

    // Get the spacepoints
    auto const spHandle = Event.getValidHandle<std::vector<recob::SpacePoint>>(fPFParticleLabel);

    // Get the hits associated with the space points
    const art::FindManyP<recob::Hit>& fmsp =
      ShowerEleHolder.GetFindManyP<recob::Hit>(spHandle, Event, fPFParticleLabel);

    //Only consider hits in the same tpcs as the vertex.
    geo::Point_t ShowerStartPosition = {-999, -999, -999};
    ShowerEleHolder.GetElement(fShowerStartPositionInputLabel, ShowerStartPosition);
    geo::TPCID vtxTPC = fGeom->FindTPCAtPosition(ShowerStartPosition);

    //Get the initial track
    recob::Track InitialTrack;
    ShowerEleHolder.GetElement(fInitialTrackInputLabel, InitialTrack);

    double pfpT0Time(0); // If no T0 found, assume the particle happened at trigger time (0)
    if (fT0Correct) {
      auto const pfpHandle = Event.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
      const art::FindManyP<anab::T0>& fmpfpt0 =
        ShowerEleHolder.GetFindManyP<anab::T0>(pfpHandle, Event, fPFParticleLabel);
      std::vector<art::Ptr<anab::T0>> pfpT0Vec = fmpfpt0.at(pfparticle.key());
      if (pfpT0Vec.size() == 1) { pfpT0Time = pfpT0Vec.front()->Time(); }
    }

    //Don't care that I could use a vector.
    std::map<int, std::vector<double>> dEdx_vec;
    std::map<int, std::vector<double>> dEdx_vecErr;
    std::map<int, int> num_hits;

    for (unsigned i = 0; i != fChannelMap.MaxPlanes(); ++i) {
      dEdx_vec[i] = {};
      dEdx_vecErr[i] = {};
      num_hits[i] = 0;
    }

    auto const clockData =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(Event);
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(Event, clockData);

    std::map<art::Ptr<recob::Hit>, std::vector<art::Ptr<recob::Hit>>> hitSnippets;
    if (fSumHitSnippets) {
      std::vector<art::Ptr<recob::Hit>> trackHits;
      ShowerEleHolder.GetElement(fInitialTrackHitsInputLabel, trackHits);

      hitSnippets = IShowerTool::GetLArPandoraShowerAlg().OrganizeHits(trackHits);
    }

    //Loop over the spacepoints
    for (auto const& sp : tracksps) {

      //Get the associated hit
      std::vector<art::Ptr<recob::Hit>> hits = fmsp.at(sp.key());
      if (hits.empty()) {
        if (fVerbose)
          mf::LogWarning("ShowerTrajPointdEdx")
            << "no hit for the spacepoint. This suggest the find many is wrong." << std::endl;
        continue;
      }
      const art::Ptr<recob::Hit> hit = hits[0];

      if (fSumHitSnippets && !hitSnippets.count(hit)) continue;

      double wirepitch = fChannelMap.Plane(hit->WireID()).WirePitch();

      //Only consider hits in the same tpc
      geo::PlaneID planeid = hit->WireID();
      geo::TPCID TPC = planeid.asTPCID();
      if (TPC != vtxTPC) { continue; }

      //Ignore spacepoints within a few wires of the vertex.
      auto const pos = sp->position();
      double dist_from_start = (pos - ShowerStartPosition).R();

      if (fCutStartPosition) {
        if (dist_from_start < fMinDistCutOff * wirepitch) { continue; }

        if (dist_from_start > dEdxTrackLength) { continue; }
      }

      //Find the closest trajectory point of the track. These should be in order if the user has used ShowerTrackTrajToSpacePoint_tool but the sake of gernicness I'll get the cloest sp.
      unsigned int index = 999;
      double MinDist = 999;
      for (unsigned int traj = 0; traj < InitialTrack.NumberTrajectoryPoints(); ++traj) {

        //ignore bogus info.
        auto flags = InitialTrack.FlagsAtPoint(traj);
        if (flags.isSet(recob::TrajectoryPointFlagTraits::NoPoint)) { continue; }

        geo::Point_t const TrajPosition = InitialTrack.LocationAtPoint(traj);

        auto const dist = (pos - TrajPosition).R();

        if (dist < MinDist && dist < MaxDist * wirepitch) {
          MinDist = dist;
          index = traj;
        }
      }

      //If there is no matching trajectory point then bail.
      if (index == 999) { continue; }

      geo::Point_t const TrajPosition = InitialTrack.LocationAtPoint(index);
      geo::Point_t const TrajPositionStart = InitialTrack.LocationAtPoint(0);

      //Ignore values with 0 mag from the start position
      if ((TrajPosition - TrajPositionStart).R() == 0) { continue; }
      if ((TrajPosition - ShowerStartPosition).R() == 0) { continue; }

      if ((TrajPosition - TrajPositionStart).R() < fMinDistCutOff * wirepitch) { continue; }

      //Get the direction of the trajectory point
      geo::Vector_t const TrajDirection = InitialTrack.DirectionAtPoint(index);

      //If the direction is in the same direction as the wires within some tolerance the hit finding struggles. Let remove these.
      // Note that we project in the YZ plane to make sure we are not cutting on
      // the angle into the wire planes, that should be done by the shaping time cut
      geo::Vector_t const TrajDirectionYZ{0, TrajDirection.Y(), TrajDirection.Z()};
      auto const PlaneDirection = fChannelMap.Plane(planeid).GetIncreasingWireDirection();

      if (std::abs((TMath::Pi() / 2 - Angle(TrajDirectionYZ, PlaneDirection))) < fMinAngleToWire) {
        if (fVerbose) mf::LogWarning("ShowerTrajPointdEdx") << "remove from angle cut" << std::endl;
        continue;
      }

      //If the direction is too much into the wire plane then the shaping amplifer cuts the charge. Lets remove these events.
      double velocity = detProp.DriftVelocity(detProp.Efield(), detProp.Temperature());
      double distance_in_x = TrajDirection.X() * (wirepitch / TrajDirection.Dot(PlaneDirection));
      double time_taken = std::abs(distance_in_x / velocity);

      //Shaping time doesn't seem to exist in a global place so add it as a fcl.
      if (fShapingTime < time_taken) {
        if (fVerbose) mf::LogWarning("ShowerTrajPointdEdx") << "move for shaping time" << std::endl;
        continue;
      }

      if ((TrajPosition - TrajPositionStart).R() > dEdxTrackLength) { continue; }

      //Iterate the number of hits on the plane
      ++num_hits[planeid.Plane];

      //If we still exist then we can be used in the calculation. Calculate the 3D pitch
      double trackpitch = (TrajDirection * (wirepitch / TrajDirection.Dot(PlaneDirection))).R();

      if (fSCECorrectPitch) {
        trackpitch = IShowerTool::GetLArPandoraShowerAlg().SCECorrectPitch(
          trackpitch, pos, TrajDirection.Unit(), hit->WireID().TPC);
      }

      //Calculate the dQdx
      double dQdx = hit->Integral();
      if (fSumHitSnippets) {
        for (const art::Ptr<recob::Hit> secondaryHit : hitSnippets[hit])
          dQdx += secondaryHit->Integral();
      }
      dQdx /= trackpitch;

      //Calculate the dEdx
      double localEField = detProp.Efield();
      if (fSCECorrectEField) {
        localEField = IShowerTool::GetLArPandoraShowerAlg().SCECorrectEField(localEField, pos);
      }
      double dEdx = fCalorimetryAlg.dEdx_AREA(
        clockData, detProp, dQdx, hit->PeakTime(), planeid.Plane, pfpT0Time, localEField);

      //Add the value to the dEdx
      dEdx_vec[planeid.Plane].push_back(dEdx);
    }

    //Choose max hits based on hitnum
    int max_hits = 0;
    int best_plane = -std::numeric_limits<int>::max();
    for (auto const& [plane, numHits] : num_hits) {
      if (fVerbose > 2) std::cout << "Plane: " << plane << " with size: " << numHits << std::endl;
      if (numHits > max_hits) {
        best_plane = plane;
        max_hits = numHits;
      }
    }

    if (best_plane < 0) {
      if (fVerbose)
        mf::LogError("ShowerTrajPointdEdx") << "No hits in any plane, returning " << std::endl;
      return 1;
    }

    //Search for blow ups and gradient changes.
    //Electrons have a very flat dEdx as function of energy till ~10MeV.
    //If there is a sudden jump particle has probably split
    //If there is very large dEdx we have either calculated it wrong (probably) or the Electron is coming to end.
    //Assumes hits are ordered!
    std::map<int, std::vector<double>> dEdx_vec_cut;
    for (geo::PlaneID plane_id : fChannelMap.Iterate<geo::PlaneID>()) {
      dEdx_vec_cut[plane_id.Plane] = {};
    }

    for (auto& dEdx_plane : dEdx_vec) {
      FinddEdxLength(dEdx_plane.second, dEdx_vec_cut[dEdx_plane.first]);
    }

    //Never have the stats to do a landau fit and get the most probable value. User decides if they want the median value or the mean.
    std::vector<double> dEdx_val;
    std::vector<double> dEdx_valErr;
    for (auto const& dEdx_plane : dEdx_vec_cut) {

      if ((dEdx_plane.second).empty()) {
        dEdx_val.push_back(-999);
        dEdx_valErr.push_back(-999);
        continue;
      }

      if (fUseMedian) {
        dEdx_val.push_back(TMath::Median((dEdx_plane.second).size(), &(dEdx_plane.second)[0]));
      }
      else {
        //Else calculate the mean value.
        double dEdx_mean = 0;
        for (auto const& dEdx : dEdx_plane.second) {
          if (dEdx > 10 || dEdx < 0) { continue; }
          dEdx_mean += dEdx;
        }
        dEdx_val.push_back(dEdx_mean / (float)(dEdx_plane.second).size());
      }
    }

    if (fVerbose > 1) {
      std::cout << "#Best Plane: " << best_plane << std::endl;
      for (unsigned int plane = 0; plane < dEdx_vec.size(); plane++) {
        std::cout << "#Plane: " << plane << " #" << std::endl;
        std::cout << "#Median: " << dEdx_val[plane] << " #" << std::endl;
        if (fVerbose > 2) {
          for (auto const& dEdx : dEdx_vec_cut[plane]) {
            std::cout << "dEdx: " << dEdx << std::endl;
          }
        }
      }
    }

    //Need to sort out errors sensibly.
    ShowerEleHolder.SetElement(dEdx_val, dEdx_valErr, fShowerdEdxOutputLabel);
    ShowerEleHolder.SetElement(best_plane, fShowerBestPlaneOutputLabel);
    ShowerEleHolder.SetElement(dEdx_vec_cut, fShowerdEdxVecOutputLabel);
    return 0;
  }

  void ShowerTrajPointdEdx::FinddEdxLength(std::vector<double>& dEdx_vec,
                                           std::vector<double>& dEdx_val)
  {

    //As default do not apply this cut.
    if (fdEdxCut > 10) {
      dEdx_val = dEdx_vec;
      return;
    }

    //Can only do this with 4 hits.
    if (dEdx_vec.size() < 4) {
      dEdx_val = dEdx_vec;
      return;
    }

    bool upperbound = false;

    //See if we are in the upper bound or upper bound defined by the cut.
    int upperbound_int = 0;
    if (dEdx_vec[0] > fdEdxCut) { ++upperbound_int; }
    if (dEdx_vec[1] > fdEdxCut) { ++upperbound_int; }
    if (dEdx_vec[2] > fdEdxCut) { ++upperbound_int; }
    if (upperbound_int > 1) { upperbound = true; }

    dEdx_val.push_back(dEdx_vec[0]);
    dEdx_val.push_back(dEdx_vec[1]);
    dEdx_val.push_back(dEdx_vec[2]);

    for (unsigned int dEdx_iter = 2; dEdx_iter < dEdx_vec.size(); ++dEdx_iter) {

      //The Function of dEdx as a function of E is flat above ~10 MeV.
      //We are looking for a jump up (or down) above the ladau width in the dEx
      //to account account for pair production.
      //Dom Estimates that the somwhere above 0.28 MeV will be a good cut but 999 will prevent this stage.
      double dEdx = dEdx_vec[dEdx_iter];

      //We are really poo at physics and so attempt to find the pair production
      if (upperbound) {
        if (dEdx > fdEdxCut) {
          dEdx_val.push_back(dEdx);
          if (fVerbose > 1) std::cout << "Adding dEdx: " << dEdx << std::endl;
          continue;
        }
        else {
          //Maybe its a landau fluctation lets try again.
          if (dEdx_iter < dEdx_vec.size() - 1) {
            if (dEdx_vec[dEdx_iter + 1] > fdEdxCut) {
              if (fVerbose > 1)
                std::cout << "Next dEdx hit is good removing hit" << dEdx << std::endl;
              continue;
            }
          }
          //I'll let one more value
          if (dEdx_iter < dEdx_vec.size() - 2) {
            if (dEdx_vec[dEdx_iter + 2] > fdEdxCut) {
              if (fVerbose > 1)
                std::cout << "Next Next dEdx hit is good removing hit" << dEdx << std::endl;
              continue;
            }
          }
          //We are hopefully we have one of our electrons has died.
          break;
        }
      }
      else {
        if (dEdx < fdEdxCut) {
          dEdx_val.push_back(dEdx);
          if (fVerbose > 1) std::cout << "Adding dEdx: " << dEdx << std::endl;
          continue;
        }
        else {
          //Maybe its a landau fluctation lets try again.
          if (dEdx_iter < dEdx_vec.size() - 1) {
            if (dEdx_vec[dEdx_iter + 1] > fdEdxCut) {
              if (fVerbose > 1)
                std::cout << "Next dEdx hit is good removing hit " << dEdx << std::endl;
              continue;
            }
          }
          //I'll let one more value
          if (dEdx_iter < dEdx_vec.size() - 2) {
            if (dEdx_vec[dEdx_iter + 2] > fdEdxCut) {
              if (fVerbose > 1)
                std::cout << "Next Next dEdx hit is good removing hit " << dEdx << std::endl;
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
