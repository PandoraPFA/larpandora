//############################################################################
//### Name:        ShowerDecisionBestPlane_tool.cc                         ###
//### Author:      Ed Tyley                                                ###
//### Date:        13.05.19                                                ###
//### Description: Decision tool to choose best plane. Choose plane with   ###
//###              the most track hits, if it is set, otherwise plane with ###
//###              most hits                                               ###
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
#include "lardataobj/RecoBase/PFParticle.h"

namespace ShowerRecoTools {


  class ShowerDecisionBestPlane: public IShowerTool {

    public:

      ShowerDecisionBestPlane(const fhicl::ParameterSet& pset);

      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder
          ) override;

    private:

      std::string fInitialTrackBestPlaneInputLabel;
      std::string fOverallBestPlaneInputLabel;
      std::string fBestPlaneOutputLabel;
  };


  ShowerDecisionBestPlane::ShowerDecisionBestPlane(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fInitialTrackBestPlaneInputLabel(pset.get<std::string>("InitialTrackBestPlaneInputLabel")),
    fOverallBestPlaneInputLabel(pset.get<std::string>("OverallBestPlaneInputLabel")),
    fBestPlaneOutputLabel(pset.get<std::string>("BestPlaneOutputLabel"))
  {
  }

  int ShowerDecisionBestPlane::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder){

    int  bestPlane = -999;

    // If the dEdx tool has decided on a best plane, use that
    // Otherwise use the plane with the most hits, done by the linear energy tool
    if (ShowerEleHolder.CheckElement(fInitialTrackBestPlaneInputLabel)){
      ShowerEleHolder.GetElement(fInitialTrackBestPlaneInputLabel, bestPlane);
      ShowerEleHolder.SetElement(bestPlane, fBestPlaneOutputLabel);
      return 0;
    } else if (ShowerEleHolder.CheckElement(fOverallBestPlaneInputLabel)){
      ShowerEleHolder.GetElement(fOverallBestPlaneInputLabel, bestPlane);
      ShowerEleHolder.SetElement(bestPlane, fBestPlaneOutputLabel);
      return 0;
    }

    // If neither are set, do not fill the element holder and return an error.
    // The LArPandoraModularShower module will take care of it.
    mf::LogError("ShowerDecisionBestPlane") << "Shower best plane not set"<< std::endl;
    return 1;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerDecisionBestPlane)
