//############################################################################
//### Name:        ShowerLinearEnergy                                      ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the Energy of the shower. Derived      ###
//###              from the linear energy algorithm, written for           ###
//###              the EMShower_module.cc                                  ###
//############################################################################

#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft Includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"

//C++ Includes
#include <iostream>
#include <vector>

//Root Includes

namespace ShowerRecoTools {

  class ShowerLinearEnergy:IShowerTool {

    public:

      ShowerLinearEnergy(const fhicl::ParameterSet& pset);

      ~ShowerLinearEnergy();

      //Physics Function. Calculate the shower Energy.
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerElementHolder
          ) override;
    private:

      double CalculateEnergy(std::vector<art::Ptr<recob::Hit> >& hits, int& plane);

      //fcl parameters
      unsigned int        fNumPlanes;
      std::vector<double> fGradients;   //Gradient of the linear fit of total charge to total energy
      std::vector<double> fIntercepts;  //Intercept of the linear fit of total charge to total energy

      art::InputTag fPFParticleModuleLabel;

      std::string fShowerEnergyOutputLabel;
      std::string fShowerBestPlaneOutputLabel;

      //Services
      detinfo::DetectorProperties const* detprop = nullptr;
      art::ServiceHandle<geo::Geometry> fGeom;

  };

  ShowerLinearEnergy::ShowerLinearEnergy(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fGradients(pset.get<std::vector<double> >("Gradients")),
    fIntercepts(pset.get<std::vector<double> >("Intercepts")),
    fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel")),
    fShowerEnergyOutputLabel(pset.get<std::string>("ShowerEnergyOutputLabel")),
    fShowerBestPlaneOutputLabel(pset.get<std::string>("ShowerBestPlaneOutputLabel")),
    detprop(lar::providerFrom<detinfo::DetectorPropertiesService>())
  {
    fNumPlanes = fGeom->Nplanes();
    if (fNumPlanes!=fGradients.size() || fNumPlanes!=fIntercepts.size()){
      throw cet::exception("ShowerLinearEnergy")
        << "The number of planes does not match the size of the fcl parametes passed: Num Planes: "
        << fNumPlanes << ", Gradients size: " << fGradients.size() << ", Intercpts size: "
        << fIntercepts.size();
    }
  }

  ShowerLinearEnergy::~ShowerLinearEnergy()
  {
  }

  int ShowerLinearEnergy::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder
      ){

    // Get the assocated pfParicle vertex PFParticles
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
      throw cet::exception("ShowerLinearEnergy") << "Could not get the pandora pf particles. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      return 1;
    }

    //Get the clusters
    art::Handle<std::vector<recob::Cluster> > clusHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, clusHandle)){
      throw cet::exception("ShowerLinearEnergy") << "Could not get the pandora clusters. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      return 1;
    }
    art::FindManyP<recob::Cluster>& fmc = ShowerEleHolder.GetFindManyP<recob::Cluster>(
        pfpHandle, Event, fPFParticleModuleLabel);
    // art::FindManyP<recob::Cluster> fmc(pfpHandle, Event, fPFParticleModuleLabel);
    std::vector<art::Ptr<recob::Cluster> > clusters = fmc.at(pfparticle.key());

    //Get the hit association
    art::FindManyP<recob::Hit>& fmhc = ShowerEleHolder.GetFindManyP<recob::Hit>(
        clusHandle, Event, fPFParticleModuleLabel);
    // art::FindManyP<recob::Hit> fmhc(clusHandle, Event, fPFParticleModuleLabel);

    std::map<unsigned int, std::vector<art::Ptr<recob::Hit> > > planeHits;

    //Loop over the clusters in the plane and get the hits
    for(auto const& cluster: clusters){

      //Get the hits
      std::vector<art::Ptr<recob::Hit> > hits = fmhc.at(cluster.key());

      //Get the plane.
      unsigned int plane = cluster->Plane().Plane;

      planeHits[plane].insert(planeHits[plane].end(),hits.begin(),hits.end());
    }

    // Calculate the energy fro each plane && best plane
    int bestPlane                 = -999;
    unsigned int bestPlaneNumHits = 0;

    //Holder for the final product
    std::vector<double> energyVec(fNumPlanes, -999);
    std::vector<double> energyError(fNumPlanes, -999);

    for(auto const& planeHitIter: planeHits){

      std::vector<art::Ptr<recob::Hit> > hits = planeHitIter.second;

      int plane                 = planeHitIter.first;
      unsigned int planeNumHits = hits.size();

      //Calculate the Energy for
      double Energy = CalculateEnergy(hits,plane);
      // If the energy is negative, leave it at -999
      if (Energy>0)
        energyVec.at(plane) = Energy;

      if (planeNumHits > bestPlaneNumHits) {
        bestPlane        = plane;
        bestPlaneNumHits = planeNumHits;
      }
    }

    //TODO

    ShowerEleHolder.SetElement(energyVec, energyError, fShowerEnergyOutputLabel);
    // Only set the best plane if it has some hits in it
    if (bestPlane!=-999){
      ShowerEleHolder.SetElement(bestPlane, fShowerBestPlaneOutputLabel);
    }

    return 0;
  }

  //Function to calculate the energy of a shower in a plane. Using a linear map between charge and Energy.
  //Exactly the same method as the ShowerEnergyAlg.cxx. Thanks Mike.
  double ShowerLinearEnergy::CalculateEnergy(std::vector<art::Ptr<recob::Hit> >& hits, int& plane) {

    double totalCharge = 0, totalEnergy = 0;

    for (auto const& hit: hits){
      totalCharge += (hit->Integral() * TMath::Exp( (detprop->SamplingRate() * hit->PeakTime()) / (detprop->ElectronLifetime()*1e3) ) );
    }

    totalEnergy = (totalCharge * fGradients.at(plane)) + fIntercepts.at(plane);

    return totalEnergy;

  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerLinearEnergy)


