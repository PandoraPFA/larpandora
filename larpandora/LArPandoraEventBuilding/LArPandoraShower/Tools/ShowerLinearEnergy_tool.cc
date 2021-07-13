//############################################################################
//### Name:        ShowerLinearEnergy                                      ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the Energy of the shower. Derived      ###
//###              from the linear energy algorithm, written for           ###
//###              the EMShower_module.cc                                  ###
//############################################################################

//Framework Includes
#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "lardataobj/RecoBase/Cluster.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"

namespace ShowerRecoTools {

  class ShowerLinearEnergy : IShowerTool {

  public:
    ShowerLinearEnergy(const fhicl::ParameterSet& pset);

    //Physics Function. Calculate the shower Energy.
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerElementHolder) override;

  private:
    double CalculateEnergy(const detinfo::DetectorClocksData& clockData,
                           const detinfo::DetectorPropertiesData& detProp,
                           const std::vector<art::Ptr<recob::Hit>>& hits,
                           const geo::PlaneID::PlaneID_t plane) const;

    //fcl parameters
    unsigned int fNumPlanes;
    std::vector<double> fGradients;  //Gradient of the linear fit of total charge to total energy
    std::vector<double> fIntercepts; //Intercept of the linear fit of total charge to total energy

    art::InputTag fPFParticleLabel;
    int fVerbose;

    std::string fShowerEnergyOutputLabel;
    std::string fShowerBestPlaneOutputLabel;

    //Services
    art::ServiceHandle<geo::Geometry> fGeom;
  };

  ShowerLinearEnergy::ShowerLinearEnergy(const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fGradients(pset.get<std::vector<double>>("Gradients"))
    , fIntercepts(pset.get<std::vector<double>>("Intercepts"))
    , fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel"))
    , fVerbose(pset.get<int>("Verbose"))
    , fShowerEnergyOutputLabel(pset.get<std::string>("ShowerEnergyOutputLabel"))
    , fShowerBestPlaneOutputLabel(pset.get<std::string>("ShowerBestPlaneOutputLabel"))
  {
    fNumPlanes = fGeom->Nplanes();
    if (fNumPlanes != fGradients.size() || fNumPlanes != fIntercepts.size()) {
      throw cet::exception("ShowerLinearEnergy")
        << "The number of planes does not match the size of the fcl parametes passed: Num Planes: "
        << fNumPlanes << ", Gradients size: " << fGradients.size()
        << ", Intercpts size: " << fIntercepts.size();
    }
  }

  int
  ShowerLinearEnergy::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                                       art::Event& Event,
                                       reco::shower::ShowerElementHolder& ShowerEleHolder)
  {

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
    std::vector<double> energyVec(fNumPlanes, -999.);
    std::vector<double> energyError(fNumPlanes, -999.);

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

  //Function to calculate the energy of a shower in a plane. Using a linear map between charge and Energy.
  //Exactly the same method as the ShowerEnergyAlg.cxx. Thanks Mike.
  double
  ShowerLinearEnergy::CalculateEnergy(const detinfo::DetectorClocksData& clockData,
                                      const detinfo::DetectorPropertiesData& detProp,
                                      const std::vector<art::Ptr<recob::Hit>>& hits,
                                      const geo::PlaneID::PlaneID_t plane) const
  {

    double totalCharge = 0, totalEnergy = 0;

    for (auto const& hit : hits) {
      totalCharge += (hit->Integral() * std::exp((sampling_rate(clockData) * hit->PeakTime()) /
                                                 (detProp.ElectronLifetime() * 1e3)));
    }

    totalEnergy = (totalCharge * fGradients.at(plane)) + fIntercepts.at(plane);

    return totalEnergy;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerLinearEnergy)
