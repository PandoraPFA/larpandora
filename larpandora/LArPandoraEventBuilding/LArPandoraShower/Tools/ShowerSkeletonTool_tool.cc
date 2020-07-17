//############################################################################
//### Name:        ShowerSkeletonTool                                      ###
//### Author:      You                                                     ###
//### Date:        13.05.19                                                ###
//### Description: Generic form of the shower tools                        ###
//############################################################################

//Framework Includes
#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"

namespace ShowerRecoTools {


  class ShowerSkeletonTool: public IShowerTool {

    public:

      ShowerSkeletonTool(const fhicl::ParameterSet& pset);

      //Generic Direction Finder
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder
          ) override;

    private:

      //Function to add the assoctions
      int AddAssociations(const art::Ptr<recob::PFParticle>& pfpPtr, art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder) override;

      // Stuff you will probably need that inherits from the module
      art::InputTag              fPFParticleLabel;
      int                        fVerbose;
  };


  ShowerSkeletonTool::ShowerSkeletonTool(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel")),
    fVerbose(pset.get<int>("Verbose"))
  {
  }

  int ShowerSkeletonTool::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder){
    return 0;
  }

  int ShowerSkeletonTool::AddAssociations(const art::Ptr<recob::PFParticle>& pfpPtr, art::Event& Event,
      reco::shower::ShowerElementHolder& ShowerEleHolder
      ){
    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerSkeletonTool)
