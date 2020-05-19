//############################################################################
//### Name:        ShowerPCAEigenvalueLength                               ###
//### Author:      Ed Tyley                                                ###
//### Date:        07.01.20                                                ###
//### Description: Simple code to calculate the lenght from the PCA        ###
//###              eigenvalues                                             ###
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
#include "lardataobj/RecoBase/PCAxis.h"

namespace ShowerRecoTools {


  class ShowerPCAEigenvalueLength: public IShowerTool {

    public:

      ShowerPCAEigenvalueLength(const fhicl::ParameterSet& pset);

      //Generic Direction Finder
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder
          ) override;

    private:

      art::InputTag fPFParticleModuleLabel;
      std::string fShowerPCAInputLabel;
      std::string fShowerLengthOutputLabel;
      std::string fShowerOpeningAngleOutputLabel;
      float fNSigma;
  };


  ShowerPCAEigenvalueLength::ShowerPCAEigenvalueLength(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel")),
    fShowerPCAInputLabel(pset.get<std::string>("ShowerPCAInputLabel")),
    fShowerLengthOutputLabel(pset.get<std::string>("ShowerLengthOutputLabel")),
    fShowerOpeningAngleOutputLabel(pset.get<std::string>("ShowerOpeningAngleOutputLabel")),
    fNSigma(pset.get<float>("NSigma"))
  {
  }

  int ShowerPCAEigenvalueLength::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder){


    if(!ShowerEleHolder.CheckElement(fShowerPCAInputLabel)){
      mf::LogError("ShowerPCAEigenvalueLength") << "PCA not set, returning "<< std::endl;
      return 1;
    }

    recob::PCAxis PCA = recob::PCAxis();
    ShowerEleHolder.GetElement(fShowerPCAInputLabel,PCA);

    const double* eigenValues = PCA.getEigenValues();

    // The PCA eigenvalues give the deviance of space points around the center
    // Take the sqrt to get std. dev and take fNSigma in each direction:
    // Call the length fNSigma x 2 x std. dev. along primary eigenvalues
    // Call the width fNSigma x 2 x std. dev. along secondary eigenvalues

    //TODO: Actually calculate the erros (Maybe fNSigma+-1Sigma?)

    //Find the length
    double primaryEigenValue = (eigenValues)[0];
    double ShowerLength = TMath::Sqrt(primaryEigenValue) * 2 * fNSigma;
    double ShowerLengthError = -9999;

    //Find the width of the shower
    double secondaryEigenValue = (eigenValues)[1];
    double ShowerWidth = TMath::Sqrt(secondaryEigenValue) * 2 * fNSigma;

    double ShowerAngle = atan(ShowerWidth/ShowerLength);
    double ShowerAngleError = -9999;

    // Fill the shower element holder
    ShowerEleHolder.SetElement(ShowerLength, ShowerLengthError, fShowerLengthOutputLabel);
    ShowerEleHolder.SetElement(ShowerAngle, ShowerAngleError, fShowerOpeningAngleOutputLabel);

    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerPCAEigenvalueLength)
