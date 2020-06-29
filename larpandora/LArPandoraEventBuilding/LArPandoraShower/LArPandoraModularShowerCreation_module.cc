//################################################################################
//### Name: LArPandoraModularShowerCreation                                    ###
//### Author: Dominic Barker and Ed Tyley (e.tyley@sheffield.ac.uk)            ###
//### Date: 15.05.19                                                           ###
//### Description: Generic Shower Charaterisation module which allows the      ###
//###              the user choose which tool to calculate shower metrics.     ###
//###              For a complete shower the tools must define (use the exact  ###
//###              name) the following  metrics in the shower property holder: ###
//###              ShowerStartPosition                                         ###
//###              ShowerDirection                                             ###
//###              ShowerEnergy                                                ###
//###              ShowerdEdx                                                  ###
//###              ShowerLength                                                ###
//###              ShowerOpeningAngle                                          ###
//################################################################################

//Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft includes
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Algs/ShowerElementHolder.hh"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Algs/ShowerProduedPtrsHolder.hh"

//Root Includes
#include "TVector3.h"

//C++ Includes
#include <vector>

namespace reco {
  namespace shower {
    class LArPandoraModularShowerCreation;
  }
}

//Class

class reco::shower::LArPandoraModularShowerCreation: public art::EDProducer {
  public:

    LArPandoraModularShowerCreation(fhicl::ParameterSet const& pset);

  private:

    void produce(art::Event& evt);

    //This function returns the art::Ptr to the data object InstanceName.
    //In the background it uses the PtrMaker which requires the element index of
    //the unique ptr (iter).
    template <class T >
      art::Ptr<T> GetProducedElementPtr(std::string InstanceName, reco::shower::ShowerElementHolder& ShowerEleHolder, int iter=-1);


    //fcl object names
    art::InputTag fPFParticleLabel;
    bool          fAllowPartialShowers;
    int           fVerbose;
    bool          fUseAllParticles;

    //tool tags which calculate the characteristics of the shower
    std::string fShowerStartPositionLabel;
    std::string fShowerDirectionLabel;
    std::string fShowerEnergyLabel;
    std::string fShowerLengthLabel;
    std::string fShowerOpeningAngleLabel;
    std::string fShowerdEdxLabel;
    std::string fShowerBestPlaneLabel;

    //fcl tools
    std::vector<std::unique_ptr<ShowerRecoTools::IShowerTool> > fShowerTools;
    std::vector<std::string>                                    fShowerToolNames;

    //map to the unique ptrs to
    reco::shower::ShowerProduedPtrsHolder uniqueproducerPtrs;
};

//This function returns the art::Ptr to the data object InstanceName.
//In the background it uses the PtrMaker which requires the element index of
//the unique ptr (iter).
template <class T >
art::Ptr<T> reco::shower::LArPandoraModularShowerCreation::GetProducedElementPtr(
    std::string InstanceName, reco::shower::ShowerElementHolder& ShowerEleHolder, int iter){

  bool check_element = ShowerEleHolder.CheckElement(InstanceName);
  if(!check_element){
    throw cet::exception("LArPandoraModularShowerCreation")
      << "To get a element that does not exist" << std::endl;
    return art::Ptr<T>();
  }

  bool check_ptr = uniqueproducerPtrs.CheckUniqueProduerPtr(InstanceName);
  if(!check_ptr){
    throw cet::exception("LArPandoraModularShowerCreation")
      << "Tried to get a ptr that does not exist" << std::endl;
    return art::Ptr<T>();
  }


  //Get the number of the shower we are on.
  int index;
  if(iter != -1){
    index = iter;
  }
  else{
    index = ShowerEleHolder.GetShowerNumber();
  }

  //Make the ptr
  art::Ptr<T> artptr = uniqueproducerPtrs.GetArtPtr<T>(InstanceName,index);
  return artptr;
}

reco::shower::LArPandoraModularShowerCreation::LArPandoraModularShowerCreation(fhicl::ParameterSet const& pset) :
  EDProducer{pset},
  fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel")),
  fAllowPartialShowers(pset.get<bool>("AllowPartialShowers")),
  fVerbose(pset.get<int>("Verbose", 0)),
  fUseAllParticles(pset.get<bool>("UseAllParticles", false)),
  fShowerStartPositionLabel(pset.get<std::string>("ShowerStartPositionLabel")),
  fShowerDirectionLabel(pset.get<std::string>("ShowerDirectionLabel")),
  fShowerEnergyLabel(pset.get<std::string>("ShowerEnergyLabel")),
  fShowerLengthLabel(pset.get<std::string>("ShowerLengthLabel")),
  fShowerOpeningAngleLabel(pset.get<std::string>("ShowerOpeningAngleLabel")),
  fShowerdEdxLabel(pset.get<std::string>("ShowerdEdxLabel")),
  fShowerBestPlaneLabel(pset.get<std::string>("ShowerBestPlaneLabel"))
{
  //Intialise the tools
  auto tool_psets = pset.get<std::vector<fhicl::ParameterSet>>("ShowerFinderTools");
  for (auto& tool_pset : tool_psets) {

    std::string tool_name = tool_pset.get<std::string>("tool_type");
    // If the PFPaticle label is not set for the tool, make it use the one for the module
    // Note we also need to set the label in the Alg, via the base tool
    if (!tool_pset.has_key("PFParticleLabel")){

      // I cannot pass an art::InputTag as it is mangled, so lets make a string instead
      std::string PFParticleLabelString(fPFParticleLabel.label()+":"+fPFParticleLabel.instance()
          +":"+fPFParticleLabel.process());

      tool_pset.put<std::string>("PFParticleLabel", PFParticleLabelString);
      fhicl::ParameterSet base_pset = tool_pset.get<fhicl::ParameterSet>("BaseTools");
      fhicl::ParameterSet alg_pset = base_pset.get<fhicl::ParameterSet>("LArPandoraShowerAlg");
      alg_pset.put<std::string>("PFParticleLabel", PFParticleLabelString);
      base_pset.put_or_replace<fhicl::ParameterSet>("LArPandoraShowerAlg", alg_pset);
      tool_pset.put_or_replace<fhicl::ParameterSet>("BaseTools", base_pset);
    }

    // If we have not explicitly set verboseness for a given tool, use global level
    if (!tool_pset.has_key("Verbose")){
      tool_pset.put<int>("Verbose", fVerbose);
    }

    fShowerTools.push_back(art::make_tool<ShowerRecoTools::IShowerTool>(tool_pset));
    fShowerToolNames.push_back(tool_name);
  }

  //  Initialise the EDProducer ptr in the tools
  std::vector<std::string> SetupTools;
  for(unsigned int i=0; i<fShowerTools.size(); ++i){
    if(std::find(SetupTools.begin(), SetupTools.end(), fShowerToolNames[i]) != SetupTools.end()){continue;}
    fShowerTools[i]->SetPtr(&producesCollector());
    fShowerTools[i]->InitaliseProducerPtr(uniqueproducerPtrs);
    fShowerTools[i]->InitialiseProducers();
  }

  //Initialise the other paramters.

  produces<std::vector<recob::Shower> >();
  produces<art::Assns<recob::Shower, recob::Hit> >();
  produces<art::Assns<recob::Shower, recob::Cluster> >();
  produces<art::Assns<recob::Shower, recob::SpacePoint> >();
  produces<art::Assns<recob::Shower, recob::PFParticle> >();

  // Output -- showers and associations with hits and clusters
  uniqueproducerPtrs.SetShowerUniqueProduerPtr(type<std::vector<recob::Shower> >(),"shower");
  uniqueproducerPtrs.SetShowerUniqueProduerPtr(type<art::Assns<recob::Shower, recob::Cluster> >(),"clusterAssociationsbase");
  uniqueproducerPtrs.SetShowerUniqueProduerPtr(type<art::Assns<recob::Shower, recob::Hit> >(),"hitAssociationsbase");
  uniqueproducerPtrs.SetShowerUniqueProduerPtr(type<art::Assns<recob::Shower, recob::SpacePoint > >(),"spShowerAssociationsbase");
  uniqueproducerPtrs.SetShowerUniqueProduerPtr(type<art::Assns<recob::Shower, recob::PFParticle> >(),"pfShowerAssociationsbase");

  uniqueproducerPtrs.PrintPtrs();
}

void reco::shower::LArPandoraModularShowerCreation::produce(art::Event& evt) {

  //Ptr makers for the products
  uniqueproducerPtrs.SetPtrMakers(evt);
  reco::shower::ShowerElementHolder showerEleHolder;

  //Get the PFParticles
  art::Handle<std::vector<recob::PFParticle> > pfpHandle;
  std::vector<art::Ptr<recob::PFParticle> > pfps;
  if (evt.getByLabel(fPFParticleLabel, pfpHandle)){
    art::fill_ptr_vector(pfps, pfpHandle);
  }
  else {
    throw cet::exception("LArPandoraModularShowerCreation")
      << "pfps not loaded. Maybe you got the module label wrong?" << std::endl;
  }

  //Handle to access the pandora hits assans
  art::Handle<std::vector<recob::Cluster> > clusterHandle;
  if (!evt.getByLabel(fPFParticleLabel,clusterHandle)){
    throw cet::exception("LArPandoraModularShowerCreation")
      << "pfp clusters are not loaded." << std::endl;
  }

  //Get the assoications to hits, clusters and spacespoints
  art::FindManyP<recob::Hit> fmh = showerEleHolder.GetFindManyP<recob::Hit>(
      clusterHandle, evt, fPFParticleLabel);
  art::FindManyP<recob::Cluster> fmcp = showerEleHolder.GetFindManyP<recob::Cluster>(
      pfpHandle, evt, fPFParticleLabel);
  art::FindManyP<recob::SpacePoint> fmspp = showerEleHolder.GetFindManyP<recob::SpacePoint>(
      pfpHandle, evt, fPFParticleLabel);
  // art::FindManyP<recob::Hit> fmh(clusterHandle, evt, fPFParticleLabel);
  // art::FindManyP<recob::Cluster> fmcp(pfpHandle, evt, fPFParticleLabel);
  // art::FindManyP<recob::SpacePoint> fmspp(pfpHandle, evt, fPFParticleLabel);

  if(!fmcp.isValid()){
    throw cet::exception("LArPandoraModularShowerCreation")
      << "Find many clusters is not valid." << std::endl;
  }
  if(!fmh.isValid()){
    throw cet::exception("LArPandoraModularShowerCreation")
      << "Find many hits is not valid." << std::endl;
  }
  if(!fmspp.isValid()){
    throw cet::exception("LArPandoraModularShowerCreation")
      << "Find many spacepoints is not valid." << std::endl;
  }

  //Holder to pass to the functions, contains the 6 properties of the shower
  // - Start Poistion
  // - Direction
  // - Initial Track
  // - Initial Track Hits
  // - Energy
  // - dEdx
  // - Length
  // - Opening Angle

  int shower_iter = 0;
  //Loop of the pf particles
  for(auto const& pfp: pfps){

    //Update the shower iterator
    showerEleHolder.SetShowerNumber(shower_iter);

    //loop only over showers unless otherwise specified
    if (!fUseAllParticles && pfp->PdgCode() != 11 && pfp->PdgCode() != 22)
      continue;

    //Get the associated hits,clusters and spacepoints
    std::vector<art::Ptr<recob::Cluster> >    showerClusters    = fmcp.at(pfp.key());
    std::vector<art::Ptr<recob::SpacePoint> > showerSpacePoints = fmspp.at(pfp.key());

    // Check the pfp has at least 1 cluster (i.e. not a pfp neutrino)
    if (!showerClusters.size())
      continue;

    //Calculate the shower properties
    //Loop over the shower tools
    int err = 0;
    unsigned int i=0;
    for(auto const& fShowerTool: fShowerTools){

      //Calculate the metric
      std::string evd_disp_append = fShowerToolNames[i]+"_iteration"+std::to_string(0) + "_" +
        this->moduleDescription().moduleLabel();

      err = fShowerTool->RunShowerTool(pfp,evt,showerEleHolder,evd_disp_append);

      if(err && fVerbose){
        mf::LogError("LArPandoraModularShowerCreation") << "Error in shower tool: " << fShowerToolNames[i]
          << " with code: " << err << std::endl;
      }
      ++i;
    }

    //If we are are not allowing partial shower check all of the things
    if(!fAllowPartialShowers){
      // If we recieved an error call from a tool return;

      // Check everything we need is in the shower element holder
      if(!showerEleHolder.CheckElement(fShowerStartPositionLabel)){
        if (fVerbose)
          mf::LogError("LArPandoraModularShowerCreation")
            << "The start position is not set in the element holder. bailing" << std::endl;
        continue;
      }
      if(!showerEleHolder.CheckElement(fShowerDirectionLabel)){
        if (fVerbose)
          mf::LogError("LArPandoraModularShowerCreation")
            << "The direction is not set in the element holder. bailing" << std::endl;
        continue;
      }
      if(!showerEleHolder.CheckElement(fShowerEnergyLabel)){
        if (fVerbose)
          mf::LogError("LArPandoraModularShowerCreation")
            << "The energy is not set in the element holder. bailing" << std::endl;
        continue;
      }
      if(!showerEleHolder.CheckElement(fShowerdEdxLabel)){
        if (fVerbose)
          mf::LogError("LArPandoraModularShowerCreation")
            << "The dEdx is not set in the element holder. bailing" << std::endl;
        continue;
      }
      if(!showerEleHolder.CheckElement(fShowerBestPlaneLabel)){
        if (fVerbose)
          mf::LogError("LArPandoraModularShowerCreation")
            << "The BestPlane is not set in the element holder. bailing" << std::endl;
        continue;
      }
      if(!showerEleHolder.CheckElement(fShowerLengthLabel)){
        if (fVerbose)
          mf::LogError("LArPandoraModularShowerCreation")
            << "The length is not set in the element holder. bailing" << std::endl;
        continue;
      }
      if(!showerEleHolder.CheckElement(fShowerOpeningAngleLabel)){
        if (fVerbose)
          mf::LogError("LArPandoraModularShowerCreation")
            << "The opening angle is not set in the element holder. bailing" << std::endl;
        continue;
      }

      //Check All of the products that have been asked to be checked.
      bool elements_are_set = showerEleHolder.CheckAllElementTags();
      if(!elements_are_set){
        if (fVerbose)
          mf::LogError("LArPandoraModularShowerCreation")
            << "Not all the elements in the property holder which should be set are not. Bailing. " << std::endl;
        continue;
      }

      ///Check all the producers
      bool producers_are_set = uniqueproducerPtrs.CheckAllProducedElements(showerEleHolder);
      if(!producers_are_set){
        if (fVerbose)
          mf::LogError("LArPandoraModularShowerCreation")
            << "Not all the elements in the property holder which are produced are not set. Bailing. " << std::endl;
        continue;
      }
    }

    //Get the properties
    TVector3            ShowerStartPosition     = {-999,-999,-999};
    TVector3            ShowerDirection         = {-999,-999,-999};
    std::vector<double> ShowerEnergy            = {-999,-999,-999};
    std::vector<double> ShowerdEdx              = {-999,-999,-999};
    int                 BestPlane               = -999;
    double              ShowerLength            = -999;
    double              ShowerOpeningAngle      = -999;

    TVector3            ShowerStartPositionErr  = {-999,-999,-999};
    TVector3            ShowerDirectionErr      = {-999,-999,-999};
    std::vector<double> ShowerEnergyErr         = {-999,-999,-999};
    std::vector<double> ShowerdEdxErr           = {-999,-999,-999};

    err = 0;
    if(showerEleHolder.CheckElement(fShowerStartPositionLabel))
      err += showerEleHolder.GetElementAndError(fShowerStartPositionLabel,ShowerStartPosition,ShowerStartPositionErr);
    if(showerEleHolder.CheckElement(fShowerDirectionLabel))
      err += showerEleHolder.GetElementAndError(fShowerDirectionLabel,ShowerDirection,ShowerDirectionErr);
    if(showerEleHolder.CheckElement(fShowerEnergyLabel))
      err += showerEleHolder.GetElementAndError(fShowerEnergyLabel,ShowerEnergy,ShowerEnergyErr);
    if(showerEleHolder.CheckElement(fShowerdEdxLabel))
      err += showerEleHolder.GetElementAndError(fShowerdEdxLabel,ShowerdEdx,ShowerdEdxErr  );
    if(showerEleHolder.CheckElement(fShowerBestPlaneLabel))
      err += showerEleHolder.GetElement(fShowerBestPlaneLabel,BestPlane);
    if(showerEleHolder.CheckElement(fShowerLengthLabel))
      err += showerEleHolder.GetElement(fShowerLengthLabel,ShowerLength);
    if(showerEleHolder.CheckElement(fShowerOpeningAngleLabel))
      err += showerEleHolder.GetElement(fShowerOpeningAngleLabel,ShowerOpeningAngle);

    if(err){
      throw cet::exception("LArPandoraModularShowerCreation")
        << "Error in LArPandoraModularShowerCreation Module. A Check on a shower property failed " << std::endl;
    }

    if(fVerbose>1){
      //Check the shower
      std::cout << "Shower Vertex: X:" << ShowerStartPosition.X() << " Y: "
        << ShowerStartPosition.Y() << " Z: " << ShowerStartPosition.Z() << std::endl;
      std::cout << "Shower Direction: X:" << ShowerDirection.X() << " Y: "
        << ShowerDirection.Y() << " Z: " << ShowerDirection.Z() << std::endl;
      std::cout << "Shower dEdx: Plane 0: " << ShowerdEdx.at(0) << " Plane 1: "
        << ShowerdEdx.at(1) << " Plane 2: " << ShowerdEdx.at(2) << std::endl;
      std::cout << "Shower Energy: Plane 0: " << ShowerEnergy.at(0) << " Plane 1: "
        << ShowerEnergy.at(1) << " Plane 2: " << ShowerEnergy.at(2) << std::endl;
      std::cout << "Shower Best Plane: " << BestPlane << std::endl;
      std::cout << "Shower Length: " << ShowerLength << std::endl;
      std::cout << "Shower Opening Angle: " << ShowerOpeningAngle << std::endl;

      //Print what has been created in the shower
      showerEleHolder.PrintElements();
    }

    //Make the shower
    recob::Shower shower = recob::Shower(ShowerDirection, ShowerDirectionErr, ShowerStartPosition, ShowerDirectionErr,
        ShowerEnergy, ShowerEnergyErr, ShowerdEdx, ShowerdEdxErr, BestPlane, util::kBogusI, ShowerLength, ShowerOpeningAngle);
    showerEleHolder.SetElement(shower,"shower");
    ++shower_iter;
    art::Ptr<recob::Shower> ShowerPtr = this->GetProducedElementPtr<recob::Shower>("shower",showerEleHolder);

    //Associate the pfparticle
    uniqueproducerPtrs.AddSingle<art::Assns<recob::Shower, recob::PFParticle>>(ShowerPtr,pfp,"pfShowerAssociationsbase");

    //Add the hits for each "cluster"
    for(auto const& cluster: showerClusters){

      //Associate the clusters
      std::vector<art::Ptr<recob::Hit> > ClusterHits = fmh.at(cluster.key());
      uniqueproducerPtrs.AddSingle<art::Assns<recob::Shower, recob::Cluster>>(ShowerPtr,cluster,"clusterAssociationsbase");

      //Associate the hits
      for(auto const& hit: ClusterHits){
        uniqueproducerPtrs.AddSingle<art::Assns<recob::Shower, recob::Hit>>(ShowerPtr, hit,"hitAssociationsbase");
      }
    }

    //Associate the spacepoints
    for(auto const& sp: showerSpacePoints){
      uniqueproducerPtrs.AddSingle<art::Assns<recob::Shower, recob::SpacePoint>>(ShowerPtr,sp,"spShowerAssociationsbase");
    }

    //Loop over the tool data products and add them.
    uniqueproducerPtrs.AddDataProducts(showerEleHolder);

    //AddAssociations
    int assn_err = 0;
    for(auto const& fShowerTool: fShowerTools){
      //AddAssociations
      assn_err += fShowerTool->AddAssociations(pfp, evt,showerEleHolder);
    }
    if(!fAllowPartialShowers && assn_err > 0){
      if (fVerbose)
        mf::LogError("LArPandoraModularShowerCreation")
          << "A association failed and not allowing partial showers. The shower will not be added to the event " << std::endl;
      continue;
    }

    //Reset the showerproperty holder.
    showerEleHolder.ClearShower();
  }

  //Put everything in the event.
  uniqueproducerPtrs.MoveAllToEvent(evt);

  //Reset the ptrs to the data products
  uniqueproducerPtrs.reset();
}

DEFINE_ART_MODULE(reco::shower::LArPandoraModularShowerCreation)
