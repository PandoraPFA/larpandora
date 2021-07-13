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
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

//C++ Includes
#include <tuple>

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
                           const geo::PlaneID::PlaneID_t plane) const;

    art::InputTag fPFParticleLabel;
    int fVerbose;

    std::string fShowerEnergyOutputLabel;
    std::string fShowerBestPlaneOutputLabel;

    //Services
    art::ServiceHandle<geo::Geometry> fGeom;
    calo::CalorimetryAlg fCalorimetryAlg;

    // Declare stuff
    double fRecombinationFactor;
  };

  ShowerNumElectronsEnergy::ShowerNumElectronsEnergy(const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel"))
    , fVerbose(pset.get<int>("Verbose"))
    , fShowerEnergyOutputLabel(pset.get<std::string>("ShowerEnergyOutputLabel"))
    , fShowerBestPlaneOutputLabel(pset.get<std::string>("ShowerBestPlaneOutputLabel"))
    , fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
    , fRecombinationFactor(pset.get<double>("RecombinationFactor"))
  {}

  int
  ShowerNumElectronsEnergy::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
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

    const art::FindManyP<recob::Cluster>& fmc =
      ShowerEleHolder.GetFindManyP<recob::Cluster>(pfpHandle, Event, fPFParticleLabel);
    // art::FindManyP<recob::Cluster> fmc(pfpHandle, Event, fPFParticleLabel);
    std::vector<art::Ptr<recob::Cluster>> clusters = fmc.at(pfparticle.key());

    //Get the hit association
    const art::FindManyP<recob::Hit>& fmhc =
      ShowerEleHolder.GetFindManyP<recob::Hit>(clusHandle, Event, fPFParticleLabel);
    // art::FindManyP<recob::Hit> fmhc(clusHandle, Event, fPFParticleLabel);

    std::map<geo::PlaneID::PlaneID_t, std::vector<art::Ptr<recob::Hit>>> planeHits;

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
    std::vector<double> energyVec(fGeom->Nplanes(), -999.);
    std::vector<double> energyError(fGeom->Nplanes(), -999.);

    auto const clockData =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(Event);
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(Event, clockData);

    for (auto const& [plane, hits] : planeHits) {

      unsigned int planeNumHits = hits.size();

      //Calculate the Energy for
      double Energy = CalculateEnergy(clockData, detProp, hits, plane);
      // If the energy is negative, leave it at -999
      if (Energy > 0) energyVec.at(plane) = Energy;

      if (planeNumHits > bestPlaneNumHits) {
        bestPlane = plane;
        bestPlaneNumHits = planeNumHits;
      }
    }

    ShowerEleHolder.SetElement(energyVec, energyError, fShowerEnergyOutputLabel);
    // Only set the best plane if it has some hits in it
    if (bestPlane < fGeom->Nplanes()) {
      // Need to cast as an int for legacy default of -999
      // have to define a new variable as we pass-by-reference when filling
      int bestPlaneVal(bestPlane);
      ShowerEleHolder.SetElement(bestPlaneVal, fShowerBestPlaneOutputLabel);
    }

    return 0;
  }

  // function to calculate the reco energy
  double
  ShowerNumElectronsEnergy::CalculateEnergy(const detinfo::DetectorClocksData& clockData,
                                            const detinfo::DetectorPropertiesData& detProp,
                                            const std::vector<art::Ptr<recob::Hit>>& hits,
                                            const geo::PlaneID::PlaneID_t plane) const
  {

    double totalCharge = 0;
    double totalEnergy = 0;
    double correctedtotalCharge = 0;
    double nElectrons = 0;

    for (auto const& hit : hits) {
      totalCharge +=
        hit->Integral() *
        fCalorimetryAlg.LifetimeCorrection(
          clockData, detProp, hit->PeakTime()); // obtain charge and correct for lifetime
    }

    // correct charge due to recombination
    correctedtotalCharge = totalCharge / fRecombinationFactor;
    // calculate # of electrons and the corresponding energy
    nElectrons = fCalorimetryAlg.ElectronsFromADCArea(correctedtotalCharge, plane);
    totalEnergy = (nElectrons / util::kGeVToElectrons) * 1000; // energy in MeV
    return totalEnergy;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerNumElectronsEnergy)
