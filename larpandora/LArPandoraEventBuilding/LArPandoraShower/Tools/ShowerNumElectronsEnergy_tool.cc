//############################################################################
//### Name:        ShowerNumElectronsEnergy                                ###
//### Author:      Tom Ham                                                 ###
//### Date:        01/04/2020                                              ###
//### Description: Tool for finding the Energy of the shower by going      ###
//###              from number of hits -> number of electrons -> energy.   ###
//###              Derived from the linear energy algorithm, written for   ###
//###              the EMShower_module.cc                                  ###
//############################################################################

//Framework Includes
#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "larcore/Geometry/WireReadout.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

// For Calorimetry normalization
#include "art/Utilities/make_tool.h"
#include "larreco/Calorimetry/INormalizeCharge.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/SpacePoint.h"

//C++ Includes
#include <tuple>

using namespace lar_pandora;

namespace ShowerRecoTools {

  class ShowerNumElectronsEnergy : IShowerTool {

  public:
    ShowerNumElectronsEnergy(const fhicl::ParameterSet& pset);

    //Physics Function. Calculate the shower Energy.
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerElementHolder) override;

  private:
    double CalculateEnergy(const detinfo::DetectorClocksData& clockData,
                           const detinfo::DetectorPropertiesData& detProp,
                           const std::vector<art::Ptr<recob::Hit>>& hits,
                           const geo::PlaneID::PlaneID_t plane,
                           const art::Event& Event,  
                           const bool applyNormalization,
                           const HitsToSpacePoints& hitsToSpacePoints,
                           const geo::Vector_t& showerPCADir) const;

    // Normalization function
    double Normalize(const double dQdx,
		     const art::Event& e,
		     const recob::Hit& h,
		     const geo::Point_t& location,
		     const geo::Vector_t& direction,
		     const double t0) const;

    art::InputTag fPFParticleLabel;
    int fVerbose;

    std::string fShowerEnergyOutputLabel;
    std::string fShowerBestPlaneOutputLabel;

    std::vector< std::unique_ptr<INormalizeCharge> > fNormalizationTools;

    //Services
    geo::WireReadoutGeom const& fChannelMap = art::ServiceHandle<geo::WireReadout>()->Get();
    calo::CalorimetryAlg fCalorimetryAlg;

    // Declare stuff
    double fRecombinationFactor;
    bool fApplyCorrectionsInNorm; // Whether to instead apply calorimetry corrections in norm.

  };

  ShowerNumElectronsEnergy::ShowerNumElectronsEnergy(const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel"))
    , fVerbose(pset.get<int>("Verbose"))
    , fShowerEnergyOutputLabel(pset.get<std::string>("ShowerEnergyOutputLabel"))
    , fShowerBestPlaneOutputLabel(pset.get<std::string>("ShowerBestPlaneOutputLabel"))
    , fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
    , fRecombinationFactor(pset.get<double>("RecombinationFactor"))
    , fApplyCorrectionsInNorm(pset.get<bool>("ApplyCorrectionsInNorm"))
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

  int ShowerNumElectronsEnergy::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                                                 art::Event& Event,
                                                 reco::shower::ShowerElementHolder& ShowerEleHolder)
  {

    if (fVerbose)
      std::cout
        << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Shower Reco Energy Tool ~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        << std::endl;

    // Get the assocated pfParicle vertex PFParticles
    auto const pfpHandle = Event.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);

    //Get the clusters
    auto const clusHandle = Event.getValidHandle<std::vector<recob::Cluster>>(fPFParticleLabel);

    geo::Vector_t showerPCADir = {-999, -999, -999};
    ShowerEleHolder.GetElement("ShowerDirection", showerPCADir);

    const art::FindManyP<recob::Cluster>& fmc =
      ShowerEleHolder.GetFindManyP<recob::Cluster>(pfpHandle, Event, fPFParticleLabel);
    // art::FindManyP<recob::Cluster> fmc(pfpHandle, Event, fPFParticleLabel);
    std::vector<art::Ptr<recob::Cluster>> clusters = fmc.at(pfparticle.key());

    //Get the hit association
    const art::FindManyP<recob::Hit>& fmhc =
      ShowerEleHolder.GetFindManyP<recob::Hit>(clusHandle, Event, fPFParticleLabel);
    // art::FindManyP<recob::Hit> fmhc(clusHandle, Event, fPFParticleLabel);

    std::map<geo::PlaneID::PlaneID_t, std::vector<art::Ptr<recob::Hit>>> planeHits;


    SpacePointVector spacePointVector;
    SpacePointsToHits spacePointsToHits;
    HitsToSpacePoints hitsToSpacePoints;
    LArPandoraHelper::CollectSpacePoints(Event, fPFParticleLabel.label(), spacePointVector, spacePointsToHits, hitsToSpacePoints);

    // Setup normalization tools
    for (auto const& nt : fNormalizationTools)
      nt->setup(Event);

    //Loop over the clusters in the plane and get the hits
    for (auto const& cluster : clusters) {

      //Get the hits
      std::vector<art::Ptr<recob::Hit>> hits = fmhc.at(cluster.key());

      //Get the plane.
      const geo::PlaneID::PlaneID_t plane(cluster->Plane().Plane);

      planeHits[plane].insert(planeHits[plane].end(), hits.begin(), hits.end());
    }

    // Calculate the energy for each plane && best plane
    geo::PlaneID::PlaneID_t bestPlane = std::numeric_limits<geo::PlaneID::PlaneID_t>::max();
    unsigned int bestPlaneNumHits = 0;

    //Holder for the final product
    std::vector<double> energyVec(fChannelMap.Nplanes(), -999.);
    std::vector<double> energyError(fChannelMap.Nplanes(), -999.);

    auto const clockData =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(Event);
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(Event, clockData);

    for (auto const& [plane, hits] : planeHits) {

      unsigned int planeNumHits = hits.size();

      //Calculate the Energy for
      double Energy = CalculateEnergy(clockData, detProp, hits, plane, Event, fApplyCorrectionsInNorm, hitsToSpacePoints, showerPCADir);

      // If the energy is negative, leave it at -999
      if (Energy > 0) energyVec.at(plane) = Energy;

      if (planeNumHits > bestPlaneNumHits) {
        bestPlane = plane;
        bestPlaneNumHits = planeNumHits;
      }
    }

    ShowerEleHolder.SetElement(energyVec, energyError, fShowerEnergyOutputLabel);
    // Only set the best plane if it has some hits in it
    if (bestPlane < fChannelMap.Nplanes()) {
      // Need to cast as an int for legacy default of -999
      // have to define a new variable as we pass-by-reference when filling
      int bestPlaneVal(bestPlane);
      ShowerEleHolder.SetElement(bestPlaneVal, fShowerBestPlaneOutputLabel);
    }

    return 0;
  }

  // function to calculate the reco energy
  double ShowerNumElectronsEnergy::CalculateEnergy(const detinfo::DetectorClocksData& clockData,
                                                   const detinfo::DetectorPropertiesData& detProp,
                                                   const std::vector<art::Ptr<recob::Hit>>& hits,
                                                   const geo::PlaneID::PlaneID_t plane,
                                                   const art::Event& Event,
                                                   const bool applyNormalization = false,
                                                   const HitsToSpacePoints& hitsToSpacePoints = HitsToSpacePoints{},
                                                   const geo::Vector_t& showerPCADir = geo::Vector_t{0, 0, 0}) const
  {

    if (applyNormalization && hitsToSpacePoints.empty()) {
      if (fVerbose) {
          mf::LogError("ShowerNumElectronsEnergy") << "No hits to space points mapping provided, returning " << std::endl;
      }
      return 1;
    }

    double totalCharge = 0;
    double totalEnergy = 0;
    double correctedtotalCharge = 0;
    double nElectrons = 0;
    double totalChargePos = 0;
    geo::Point_t chargeWeightedPosition = {0, 0, 0}; // Initialize charge weighted position

    for (auto const& hit : hits) {
      totalCharge +=
        hit->Integral() *
        fCalorimetryAlg.LifetimeCorrection(
          clockData, detProp, hit->PeakTime()); // obtain charge and correct for lifetime
        
          if ( applyNormalization ) {
            HitsToSpacePoints::const_iterator hIter = hitsToSpacePoints.find(hit);
            if (hitsToSpacePoints.end() != hIter){
              const art::Ptr<recob::SpacePoint> spacepoint = hIter->second;            
              auto const& pos = spacepoint->position();  // this is a geo::Point_t
              chargeWeightedPosition += geo::Vector_t{pos.X(), pos.Y(), pos.Z()} * hit->Integral();
              totalChargePos += hit->Integral(); // Accumulate total charge
            }
        }
    }

    // correct charge due to recombination
    correctedtotalCharge = totalCharge / fRecombinationFactor;

    //std::cout << "Applying normalization: " << applyNormalization << std::endl;

    // apply normalization if needed
    if ( applyNormalization && hits.size() > 0 ) {
      //std::cout << "\t" << "Total charge before norm: " << correctedtotalCharge << std::endl;
      if (totalChargePos > 0) chargeWeightedPosition /= totalChargePos; // Normalize by total charge
      correctedtotalCharge = Normalize( correctedtotalCharge,
        Event,
        *hits.at(0),
        chargeWeightedPosition,
        showerPCADir,
        0 );
      //std::cout << "\t" << "Total charge after norm: " << correctedtotalCharge << std::endl;
    }
    // calculate # of electrons and the corresponding energy
    nElectrons = fCalorimetryAlg.ElectronsFromADCArea(correctedtotalCharge, plane);
    totalEnergy = (nElectrons / util::kGeVToElectrons) * 1000; // energy in MeV
    return totalEnergy;
  }

  double ShowerNumElectronsEnergy::Normalize(const double dQdx,
					const art::Event& e,
					const recob::Hit& h,
					const geo::Point_t& location,
					const geo::Vector_t& direction,
					const double t0) const
  {
    double ret = dQdx;
    for (auto const& nt : fNormalizationTools) {
      ret = nt->Normalize(ret, e, h, location, direction, t0);
      //std::cout << "\t norm: dQdx = " << ret << std::endl;
    }
    
    return ret;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerNumElectronsEnergy)
