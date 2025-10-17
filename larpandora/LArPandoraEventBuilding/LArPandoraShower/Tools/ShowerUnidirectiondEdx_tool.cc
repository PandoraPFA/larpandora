//############################################################################
//### Name:        ShowerUnidirectiondEdx                                  ###
//### Author:      Ed Tyley                                                ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the dEdx of the start track of the     ###
//###              shower using the standard calomitry module. Derived     ###
//###              from the EMShower_module.cc                             ###
//############################################################################

//Framework Includes
#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

// For Calorimetry normalization
#include "art/Utilities/make_tool.h"
#include "larreco/Calorimetry/INormalizeCharge.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"

using namespace lar_pandora;

namespace ShowerRecoTools {

  class ShowerUnidirectiondEdx : IShowerTool {

  public:
    ShowerUnidirectiondEdx(const fhicl::ParameterSet& pset);

    //Generic Direction Finder
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerEleHolder) override;

  private:

    // Normalization function
    const double Normalize(const double dQdx,
		     const art::Event& e,
		     const recob::Hit& h,
		     const geo::Point_t& location,
		     const geo::Vector_t& direction,
		     const double t0);

    //Define the services and algorithms
    art::ServiceHandle<geo::Geometry> fGeom;
    geo::WireReadoutGeom const& fChannelMap = art::ServiceHandle<geo::WireReadout>()->Get();
    calo::CalorimetryAlg fCalorimetryAlg;

    std::vector< std::unique_ptr<INormalizeCharge> > fNormalizationTools;

    //fcl parameters.
    int fVerbose;
    double fdEdxTrackLength,
      dEdxTrackLength;    //Max length from a hit can be to the start point in cm.
    bool fMaxHitPlane;    //Set the best planes as the one with the most hits
    bool fMissFirstPoint; //Do not use any hits from the first wire.
    bool fSumHitSnippets; // Whether to treat hits individually or only one hit per snippet
    bool fApplyCorrectionsInNorm; // Whether to instead apply calorimetry corrections in norm.

    std::string fShowerStartPositionInputLabel;
    std::string fInitialTrackHitsInputLabel;
    std::string fShowerDirectionInputLabel;
    std::string fShowerdEdxOutputLabel;
    std::string fShowerBestPlaneOutputLabel;
    art::InputTag fPFParticleLabel;
  };

  ShowerUnidirectiondEdx::ShowerUnidirectiondEdx(const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
    , fVerbose(pset.get<int>("Verbose"))
    , fdEdxTrackLength(pset.get<float>("dEdxTrackLength"))
    , fMaxHitPlane(pset.get<bool>("MaxHitPlane"))
    , fMissFirstPoint(pset.get<bool>("MissFirstPoint"))
    , fSumHitSnippets(pset.get<bool>("SumHitSnippets"))
    , fApplyCorrectionsInNorm(pset.get<bool>("ApplyCorrectionsInNorm"))
    , fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel"))
    , fInitialTrackHitsInputLabel(pset.get<std::string>("InitialTrackHitsInputLabel"))
    //, fInitialTrackInputLabel(pset.get<std::string>("InitialTrackInputLabel"))
    , fShowerDirectionInputLabel(pset.get<std::string>("ShowerDirectionInputLabel"))
    , fShowerdEdxOutputLabel(pset.get<std::string>("ShowerdEdxOutputLabel"))
    , fShowerBestPlaneOutputLabel(pset.get<std::string>("ShowerBestPlaneOutputLabel"))
    , fPFParticleLabel(pset.get<std::string>("PFParticleLabel"))

  {
    if ( fApplyCorrectionsInNorm ) {
      auto tool_psets = pset.get< std::vector< fhicl::ParameterSet > >("NormTools");

      int tCounter = 0;
      for ( auto const& tool_pset : tool_psets ) {
        //std::cout << "pushing back tools..." << tCounter << std::endl;
        tCounter++;
	      fNormalizationTools.push_back( art::make_tool<INormalizeCharge>(tool_pset) );
      }
    }
  }

  int ShowerUnidirectiondEdx::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                                               art::Event& Event,
                                               reco::shower::ShowerElementHolder& ShowerEleHolder)
  {

    dEdxTrackLength = fdEdxTrackLength;

    SpacePointVector spacePointVector;
    SpacePointsToHits spacePointsToHits;
    HitsToSpacePoints hitsToSpacePoints;
    LArPandoraHelper::CollectSpacePoints(Event, fPFParticleLabel.label(), spacePointVector, spacePointsToHits, hitsToSpacePoints);
    // std::cout << "There are " << hitsToSpacePoints.size() << " hits associated with space points" << std::endl;

    // Shower dEdx calculation
    if (!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)) {
      if (fVerbose)
        mf::LogError("ShowerUnidirectiondEdx") << "Start position not set, returning " << std::endl;
      return 1;
    }
    if (!ShowerEleHolder.CheckElement(fInitialTrackHitsInputLabel)) {
      if (fVerbose)
        mf::LogError("ShowerUnidirectiondEdx")
          << "Initial Track Hits not set returning" << std::endl;
      return 1;
    }
    if (!ShowerEleHolder.CheckElement(fShowerDirectionInputLabel)) {
      if (fVerbose)
        mf::LogError("ShowerUnidirectiondEdx") << "Shower Direction not set" << std::endl;
      return 1;
    }

    // auto const pfpHandle = Event.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
    // const art::FindManyP<recob::SpacePoint>& fmspp =
    //  ShowerEleHolder.GetFindManyP<recob::SpacePoint>(pfpHandle, Event, fPFParticleLabel);

    //Get the initial track hits
    std::vector<art::Ptr<recob::Hit>> trackhits;
    ShowerEleHolder.GetElement(fInitialTrackHitsInputLabel, trackhits);

    if (trackhits.empty()) {
      if (fVerbose)
        mf::LogWarning("ShowerUnidirectiondEdx") << "Not Hits in the initial track" << std::endl;
      return 0;
    }

    geo::Point_t ShowerStartPosition = {-999, -999, -999};
    ShowerEleHolder.GetElement(fShowerStartPositionInputLabel, ShowerStartPosition);

    geo::Vector_t showerDir = {-999, -999, -999};
    ShowerEleHolder.GetElement(fShowerDirectionInputLabel, showerDir);

    geo::Vector_t showerPCADir = {-999, -999, -999};
    ShowerEleHolder.GetElement("ShowerDirection", showerPCADir);

    //std::cout << "Shower direction from PCA: " << showerPCADir.X() << " " << showerPCADir.Y() << " " << showerPCADir.Z() << std::endl;

    //std::cout << TString(Form("Shower position %f %f %f, and direction %f %f %f", ShowerStartPosition.X(), ShowerStartPosition.Y(), ShowerStartPosition.Z(),
    //                      showerDir.X(), showerDir.Y(), showerDir.Z())) << std::endl;

    geo::TPCID vtxTPC = fGeom->FindTPCAtPosition(geo::vect::toPoint(ShowerStartPosition));

    // Split the track hits per plane
    std::vector<double> dEdxVec;
    std::vector<std::vector<art::Ptr<recob::Hit>>> trackHits;
    unsigned int numPlanes = fChannelMap.Nplanes();
    trackHits.resize(numPlanes);

    // Loop over the track hits and split into planes
    for (unsigned int hitIt = 0; hitIt < trackhits.size(); ++hitIt) {
      art::Ptr<recob::Hit> hit = trackhits.at(hitIt);
      geo::PlaneID hitWire = hit->WireID();
      geo::TPCID TPC = hitWire.asTPCID();

      //only get hits from the same TPC as the vertex
      if (TPC == vtxTPC) { (trackHits.at(hitWire.Plane)).push_back(hit); }
    }

    int bestHitsPlane = 0;
    int bestPlaneHits = 0;
    int bestPlane = -999;
    double minPitch = 999;

    auto const clockData =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(Event);
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(Event, clockData);

    for (unsigned int plane = 0; plane < numPlanes; ++plane) {
      std::vector<art::Ptr<recob::Hit>> trackPlaneHits = trackHits.at(plane);

      std::map<art::Ptr<recob::Hit>, std::vector<art::Ptr<recob::Hit>>> hitSnippets;
      if (fSumHitSnippets)
        hitSnippets = IShowerTool::GetLArPandoraShowerAlg().OrganizeHits(trackPlaneHits);

      if (trackPlaneHits.size()) {

        double dEdx = -999;
        double totQ = 0;
        double avgT = 0;
        double pitch = 0;

        //Calculate the pitch
        double wirepitch = fChannelMap.Plane(trackPlaneHits.at(0)->WireID()).WirePitch();
        double angleToVert =
          fChannelMap.WireAngleToVertical(fChannelMap.Plane(geo::PlaneID{0, 0, plane}).View(),
                                          trackPlaneHits[0]->WireID().planeID()) -
          0.5 * TMath::Pi();
        double cosgamma =
          std::abs(std::sin(angleToVert) * showerDir.Y() + std::cos(angleToVert) * showerDir.Z());

        pitch = wirepitch / cosgamma;

        if (pitch) { // Check the pitch is calculated correctly
          int nhits = 0;
          std::vector<float> vQ;

          //Get the first wire
          int w0 = trackPlaneHits.at(0)->WireID().Wire;

          geo::Point_t chargeWeightedPosition = {0, 0, 0}; // Initialize charge weighted position
          double totalCharge = 0; // Initialize total charge

          for (auto const& hit : trackPlaneHits) {

            if (fSumHitSnippets && !hitSnippets.count(hit)) continue;

            // Get the wire for each hit
            int w1 = hit->WireID().Wire;
            if (fMissFirstPoint && w0 == w1) { continue; }

            // Ignore hits that are too far away.
            if (std::abs((w1 - w0) * pitch) < dEdxTrackLength) {

              double q = hit->Integral();
              if (fSumHitSnippets) {
                for (const art::Ptr<recob::Hit> secondaryHit : hitSnippets[hit])
                  q += secondaryHit->Integral();
              }

              vQ.push_back(q);
              totQ += hit->Integral();
              avgT += hit->PeakTime();
              ++nhits;
              
              //std::cout << "There are " << fmspp.size() << " spacepoint collections" << std::endl;

              /*if (fmspp.size() <= hit.key()) {
                std::cout << "No spacepoints for this hit" << std::endl;
                continue;
              }
              std::vector<art::Ptr<recob::SpacePoint>> sps = fmspp.at(hit.key());
              std::cout << "Found " << sps.size() << " spacepoints for this hit." << std::endl;*/

              HitsToSpacePoints::const_iterator hIter = hitsToSpacePoints.find(hit);
              if (hitsToSpacePoints.end() != hIter){
                const art::Ptr<recob::SpacePoint> spacepoint = hIter->second;
                //const double X(spacepoint->XYZ()[0]);
                //std::cout << "Found spacepoint at X: " << X << std::endl;

                auto const& pos = spacepoint->position();  // this is a geo::Point_t
                chargeWeightedPosition += geo::Vector_t{pos.X(), pos.Y(), pos.Z()} * q;
                totalCharge += q; // Accumulate total charge
                //std::cout << "updating charge weighted position with q: " << q << std::endl;
              }

            }
          }

          // Calculate the final charge weighted average position
          if (totalCharge > 0) {
            chargeWeightedPosition /= totalCharge; // Normalize by total charge
            // std::cout << "Charge weighted position: (" << chargeWeightedPosition.X() << ", "
            //           << chargeWeightedPosition.Y() << ", " << chargeWeightedPosition.Z() << ")" << std::endl;
          }

          if (totQ) {
            // Check if the pitch is the smallest yet (best plane)
            if (pitch < minPitch) {
              minPitch = pitch;
              bestPlane = plane;
            }

            //Get the median and calculate the dEdx using the algorithm.
            double dQdx = TMath::Median(vQ.size(), &vQ[0]) / pitch;
            const auto& hit = trackPlaneHits.at(0);
            double dQdxNorm = dQdx;
            // Attempt the normalization //Mike 
            if ( fApplyCorrectionsInNorm ) {
              geo::Vector_t displacement(0, 0, 0);
              //std::cout << "Running the CorrectionsInNorm for showers" << std::endl;
              dQdxNorm = Normalize( dQdx,
                Event,
                *hit,
                chargeWeightedPosition,
                displacement,
                0 );
            }
            std::cout << "Unidirection: dQdx: " << dQdx << " dQdxNorm: " << dQdxNorm << std::endl;

            dEdx = fCalorimetryAlg.dEdx_AREA(
              clockData, detProp, dQdxNorm, avgT / nhits, trackPlaneHits.at(0)->WireID().Plane);

            if (isinf(dEdx)) { dEdx = -999; };

            if (nhits > bestPlaneHits || ((nhits == bestPlaneHits) && (pitch < minPitch))) {
              bestHitsPlane = plane;
              bestPlaneHits = nhits;
            }
          }
          dEdxVec.push_back(dEdx);
        }
        else {
          throw cet::exception("ShowerUnidirectiondEdx")
            << "pitch is 0. I can't think how it is 0? Stopping so I can tell you" << std::endl;
        }
      }
      else { // if not (trackPlaneHits.size())
        dEdxVec.push_back(-999);
      }
      trackPlaneHits.clear();
    } //end loop over planes

    //TODO
    std::vector<double> dEdxVecErr = {-999, -999, -999};

    ShowerEleHolder.SetElement(dEdxVec, dEdxVecErr, fShowerdEdxOutputLabel);

    //Set The best plane
    if (fMaxHitPlane) { bestPlane = bestHitsPlane; }

    if (bestPlane == -999) {
      throw cet::exception("ShowerUnidirectiondEdx") << "No best plane set";
    }
    else {
      ShowerEleHolder.SetElement(bestPlane, fShowerBestPlaneOutputLabel);
    }

    if (fVerbose > 1) {
      std::cout << "Best Plane: " << bestPlane << std::endl;
      for (unsigned int plane = 0; plane < dEdxVec.size(); plane++) {
        std::cout << "Plane: " << plane << " with dEdx: " << dEdxVec[plane] << std::endl;
      }
    }

    return 0;
  }

  const double ShowerUnidirectiondEdx::Normalize(const double dQdx,
					const art::Event& e,
					const recob::Hit& h,
					const geo::Point_t& location,
					const geo::Vector_t& direction,
					const double t0)
  {
    double ret = dQdx;
    for (auto const& nt : fNormalizationTools) {
      ret = nt->Normalize(ret, e, h, location, direction, t0);
      std::cout << "\t norm: dQdx = " << ret << std::endl;
    }
    
    return ret;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerUnidirectiondEdx)
