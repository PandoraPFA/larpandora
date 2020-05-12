//################################################################################
//### Name: LArPandoraModularShower                                                      ###
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
    class LArPandoraModularShower;
  }
}

//Class

class reco::shower::LArPandoraModularShower: public art::EDProducer {
  public:

    LArPandoraModularShower(fhicl::ParameterSet const& pset);

  private:

    void produce(art::Event& evt);

    //This function returns the art::Ptr to the data object InstanceName. In the background it uses the PtrMaker which requires the element index of
    //the unique ptr (iter).
    template <class T >
      art::Ptr<T> GetProducedElementPtr(std::string InstanceName, reco::shower::ShowerElementHolder& ShowerEleHolder, int iter=-1);


    //fcl object names
    art::InputTag fPFParticleModuleLabel;
    bool          fSecondInteration;
    bool          fAllowPartialShowers;
    bool          fVerbose;

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

//This function returns the art::Ptr to the data object InstanceName. In the background it uses the PtrMaker which requires the element index of
//the unique ptr (iter).
template <class T >
art::Ptr<T> reco::shower::LArPandoraModularShower::GetProducedElementPtr(std::string InstanceName, reco::shower::ShowerElementHolder& ShowerEleHolder, int iter){

  bool check_element = ShowerEleHolder.CheckElement(InstanceName);
  if(!check_element){
    throw cet::exception("LArPandoraModularShower") << "To get a element that does not exist" << std::endl;
    return art::Ptr<T>();
  }

  bool check_ptr = uniqueproducerPtrs.CheckUniqueProduerPtr(InstanceName);
  if(!check_ptr){
    throw cet::exception("LArPandoraModularShower") << "Tried to get a ptr that does not exist" << std::endl;
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



reco::shower::LArPandoraModularShower::LArPandoraModularShower(fhicl::ParameterSet const& pset) :
  EDProducer{pset}
{
  //Intialise the tools
  auto const tool_psets = pset.get<std::vector<fhicl::ParameterSet>>("ShowerFinderTools");
  for (auto const& tool_pset : tool_psets) {
    fShowerTools.push_back(art::make_tool<ShowerRecoTools::IShowerTool>(tool_pset));
    std::string tool_name = tool_pset.get<std::string>("tool_type");
    fShowerToolNames.push_back(tool_name);
    std::cout<< "Tools List: " << tool_name << std::endl;
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
  fPFParticleModuleLabel      = pset.get<art::InputTag>("PFParticleModuleLabel","pandora");
  fShowerStartPositionLabel   = pset.get<std::string  >("ShowerStartPositionLabel");
  fShowerDirectionLabel       = pset.get<std::string  >("ShowerDirectionLabel");
  fShowerEnergyLabel          = pset.get<std::string  >("ShowerEnergyLabel");
  fShowerLengthLabel          = pset.get<std::string  >("ShowerLengthLabel");
  fShowerOpeningAngleLabel    = pset.get<std::string  > ("ShowerOpeningAngleLabel");
  fShowerdEdxLabel            = pset.get<std::string  >("ShowerdEdxLabel");
  fShowerBestPlaneLabel       = pset.get<std::string  >("ShowerBestPlaneLabel");
  fSecondInteration           = pset.get<bool         >("SecondInteration",false);
  fAllowPartialShowers        = pset.get<bool         >("AllowPartialShowers",false);
  fVerbose                    = pset.get<bool         >("Verbose",false);

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

void reco::shower::LArPandoraModularShower::produce(art::Event& evt) {

  //Ptr makers for the products
  uniqueproducerPtrs.SetPtrMakers(evt);
  reco::shower::ShowerElementHolder selement_holder;

  //Get the PFParticles
  art::Handle<std::vector<recob::PFParticle> > pfpHandle;
  std::vector<art::Ptr<recob::PFParticle> > pfps;
  if (evt.getByLabel(fPFParticleModuleLabel, pfpHandle)){
    art::fill_ptr_vector(pfps, pfpHandle);
  }
  else {
    throw cet::exception("LArPandoraModularShower") << "pfps not loaded. Maybe you got the module label wrong?" << std::endl;
  }

  //Handle to access the pandora hits assans
  art::Handle<std::vector<recob::Cluster> > clusterHandle;
  if (!evt.getByLabel(fPFParticleModuleLabel,clusterHandle)){
    throw cet::exception("LArPandoraModularShower") << "pfp clusters are not loaded." << std::endl;
  }

  //Get the assoications to hits, clusters and spacespoints
  art::FindManyP<recob::Hit> fmh = selement_holder.GetFindManyP<recob::Hit>(
      clusterHandle, evt, fPFParticleModuleLabel);
  art::FindManyP<recob::Cluster> fmcp = selement_holder.GetFindManyP<recob::Cluster>(
      pfpHandle, evt, fPFParticleModuleLabel);
  art::FindManyP<recob::SpacePoint> fmspp = selement_holder.GetFindManyP<recob::SpacePoint>(
      pfpHandle, evt, fPFParticleModuleLabel);
  // art::FindManyP<recob::Hit> fmh(clusterHandle, evt, fPFParticleModuleLabel);
  // art::FindManyP<recob::Cluster> fmcp(pfpHandle, evt, fPFParticleModuleLabel);
  // art::FindManyP<recob::SpacePoint> fmspp(pfpHandle, evt, fPFParticleModuleLabel);

  if(!fmcp.isValid()){
    throw cet::exception("LArPandoraModularShower") << "Find many clusters is not valid." << std::endl;
  }
  if(!fmh.isValid()){
    throw cet::exception("LArPandoraModularShower") << "Find many hits is not valid." << std::endl;
  }
  if(!fmspp.isValid()){
    throw cet::exception("LArPandoraModularShower") << "Find many spacepoints is not valid." << std::endl;
  }

  //Holder to pass to the functions, contains the 6 properties of the shower
  // - Start Poistion
  // - Direction
  // - Initial Track
  // - Initial Track Hits
  // - Energy
  // - dEdx

  int shower_iter = 0;
  //Loop of the pf particles
  for(auto const& pfp: pfps){

    //Update the shower iterator
    selement_holder.SetShowerNumber(shower_iter);

    //loop only over showers.
    if(pfp->PdgCode() != 11 && pfp->PdgCode() != 22){continue;}

    //Calculate the shower properties
    //Loop over the shower tools
    int err = 0;
    unsigned int i=0;
    for(auto const& fShowerTool: fShowerTools){

      //Calculate the metric
      std::string evd_disp_append = fShowerToolNames[i]+"_iteration"+std::to_string(0) + "_" + this->moduleDescription().moduleLabel();

      err = fShowerTool->RunShowerTool(pfp,evt,selement_holder,evd_disp_append);

      if(err){
        mf::LogError("LArPandoraModularShower") << "Error in shower tool: " << fShowerToolNames[i]  << " with code: " << err << std::endl;
        break;
      }
      ++i;
    }
    //Should we do a second interaction now we have done a first pass of the calculation
    i=0;
    if(fSecondInteration){

      for(auto const& fShowerTool: fShowerTools){
        //Calculate the metric
        std::string evd_disp_append = fShowerToolNames[i]+"_iteration"+std::to_string(1) + "_" + this->moduleDescription().moduleLabel();

        err = fShowerTool->RunShowerTool(pfp,evt,selement_holder,evd_disp_append);

        if(err){
          mf::LogError("LArPandoraModularShower") << "Error in shower tool: " << fShowerToolNames[i]  << " with code: " << err << std::endl;
          break;
        }
        ++i;
      }
    }

    //If we are are not allowing partial shower check all of the things
    if(!fAllowPartialShowers){
      // If we recieved an error call from a tool return;
      if(err){
        mf::LogError("LArPandoraModularShower") << "Error on tool. Assuming all the shower products and properties were not set and bailing." << std::endl;
        continue;
      }

      // Check everything we need is in the shower element holder
      if(!selement_holder.CheckElement(fShowerStartPositionLabel)){
        mf::LogError("LArPandoraModularShower") << "The start position is not set in the element holder. bailing" << std::endl;
        continue;
      }
      if(!selement_holder.CheckElement(fShowerDirectionLabel)){
        mf::LogError("LArPandoraModularShower") << "The direction is not set in the element holder. bailing" << std::endl;
        continue;
      }
      if(!selement_holder.CheckElement(fShowerEnergyLabel)){
        mf::LogError("LArPandoraModularShower") << "The energy is not set in the element holder. bailing" << std::endl;
        continue;
      }
      if(!selement_holder.CheckElement(fShowerdEdxLabel)){
        mf::LogError("LArPandoraModularShower") << "The dEdx is not set in the element holder. bailing" << std::endl;
        continue;
      }
      if(!selement_holder.CheckElement(fShowerBestPlaneLabel)){
        mf::LogError("LArPandoraModularShower") << "The BestPlane is not set in the element holder. bailing" << std::endl;
        continue;
      }
      if(!selement_holder.CheckElement(fShowerLengthLabel)){
        mf::LogError("LArPandoraModularShower") << "The length is not set in the element holder. bailing" << std::endl;
        continue;
      }
      if(!selement_holder.CheckElement(fShowerOpeningAngleLabel)){
        mf::LogError("LArPandoraModularShower") << "The opening angle is not set in the element holder. bailing" << std::endl;
        continue;
      }

      //Check All of the products that have been asked to be checked.
      bool elements_are_set = selement_holder.CheckAllElementTags();
      if(!elements_are_set){
        mf::LogError("LArPandoraModularShower") << "Not all the elements in the property holder which should be set are not. Bailing. " << std::endl;
        continue;
      }

      ///Check all the producers
      bool producers_are_set = uniqueproducerPtrs.CheckAllProducedElements(selement_holder);
      if(!producers_are_set){
        mf::LogError("LArPandoraModularShower") << "Not all the elements in the property holder which are produced are not set. Bailing. " << std::endl;
        continue;
      }
    }

    //Get the properties
    TVector3                           ShowerStartPosition  = {-999,-999,-999};
    TVector3                           ShowerDirection      = {-999,-999,-999};
    std::vector<double>                ShowerEnergy         = {-999,-999,-999};
    std::vector<double>                ShowerdEdx           = {-999,-999,-999};
    int                                BestPlane            = -999;
    double                             ShowerLength         = -999;
    double                             ShowerOpeningAngle   = -999;

    TVector3                           ShowerStartPositionErr  = {-999,-999,-999};
    TVector3                           ShowerDirectionErr      = {-999,-999,-999};
    std::vector<double>                ShowerEnergyErr         = {-999,-999,-999};
    std::vector<double>                ShowerdEdxErr           = {-999,-999,-999};

    err = 0;
    if(selement_holder.CheckElement(fShowerStartPositionLabel))    err += selement_holder.GetElementAndError(fShowerStartPositionLabel,ShowerStartPosition,ShowerStartPositionErr);
    if(selement_holder.CheckElement(fShowerDirectionLabel))        err += selement_holder.GetElementAndError(fShowerDirectionLabel,ShowerDirection,ShowerDirectionErr);
    if(selement_holder.CheckElement(fShowerEnergyLabel))           err += selement_holder.GetElementAndError(fShowerEnergyLabel,ShowerEnergy,ShowerEnergyErr);
    if(selement_holder.CheckElement(fShowerdEdxLabel))             err += selement_holder.GetElementAndError(fShowerdEdxLabel,ShowerdEdx,ShowerdEdxErr  );
    if(selement_holder.CheckElement(fShowerBestPlaneLabel))        err += selement_holder.GetElement(fShowerBestPlaneLabel,BestPlane);
    if(selement_holder.CheckElement(fShowerLengthLabel))           err += selement_holder.GetElement(fShowerLengthLabel,ShowerLength);
    if(selement_holder.CheckElement(fShowerOpeningAngleLabel))     err += selement_holder.GetElement(fShowerOpeningAngleLabel,ShowerOpeningAngle);

    if(err){
      throw cet::exception("LArPandoraModularShower")  << "Error in LArPandoraModularShower Module. A Check on a shower property failed " << std::endl;
    }

    if(fVerbose){
      //Check the shower
      std::cout<<"Shower Vertex: X:"<<ShowerStartPosition.X()<<" Y: "<<ShowerStartPosition.Y()<<" Z: "<<ShowerStartPosition.Z()<<std::endl;
      std::cout<<"Shower Direction: X:"<<ShowerDirection.X()<<" Y: "<<ShowerDirection.Y()<<" Z: "<<ShowerDirection.Z()<<std::endl;
      std::cout<<"Shower dEdx: size: "<<ShowerdEdx.size()<<" Plane 0: "<<ShowerdEdx.at(0)<<" Plane 1: "<<ShowerdEdx.at(1)<<" Plane 2: "<<ShowerdEdx.at(2)<<std::endl;
      std::cout<<"Shower Energy: size: "<<ShowerEnergy.size()<<" Plane 0: "<<ShowerEnergy.at(0)<<" Plane 1: "<<ShowerEnergy.at(1)<<" Plane 2: "<<ShowerEnergy.at(2)<<std::endl;
      std::cout<<"Shower Best Plane: "<<BestPlane<<std::endl;
      std::cout<<"Shower Length: " << ShowerLength << std::endl;
      std::cout<<"Shower Opening Angle: " << ShowerOpeningAngle << std::endl;

      //Print what has been created in the shower
      selement_holder.PrintElements();
    }

    //Make the shower
    recob::Shower shower = recob::Shower(ShowerDirection, ShowerDirectionErr,ShowerStartPosition, ShowerDirectionErr,ShowerEnergy,ShowerEnergyErr,ShowerdEdx, ShowerdEdxErr, BestPlane,util::kBogusI, ShowerLength, ShowerOpeningAngle);
    selement_holder.SetElement(shower,"shower");
    ++shower_iter;
    art::Ptr<recob::Shower> ShowerPtr = this->GetProducedElementPtr<recob::Shower>("shower",selement_holder);

    //Associate the pfparticle
    uniqueproducerPtrs.AddSingle<art::Assns<recob::Shower, recob::PFParticle>>(ShowerPtr,pfp,"pfShowerAssociationsbase");

    //Get the associated hits,clusters and spacepoints
    std::vector<art::Ptr<recob::Cluster> >    showerClusters    = fmcp.at(pfp.key());
    std::vector<art::Ptr<recob::SpacePoint> > showerSpacePoints = fmspp.at(pfp.key());

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
    uniqueproducerPtrs.AddDataProducts(selement_holder);

    //AddAssociations
    int assn_err = 0;
    for(auto const& fShowerTool: fShowerTools){
      //AddAssociations
      assn_err += fShowerTool->AddAssociations(pfp, evt,selement_holder);
    }
    if(!fAllowPartialShowers && assn_err > 0){
      mf::LogError("LArPandoraModularShower") << "A association failed and you are not allowing partial showers. The event will not be added to the event " << std::endl;
      continue;
    }

    //Reset the showerproperty holder.
    selement_holder.ClearShower();
  }

  //Put everything in the event.
  uniqueproducerPtrs.MoveAllToEvent(evt);

  //Reset the ptrs to the data products
  uniqueproducerPtrs.reset();

}

DEFINE_ART_MODULE(reco::shower::LArPandoraModularShower)
