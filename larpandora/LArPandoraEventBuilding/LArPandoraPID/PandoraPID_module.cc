////////////////////////////////////////////////////////////////////////
// Class:       PandoraPID
// Plugin Type: producer (Unknown Unknown)
// File:        PandoraPID_module.cc
//
// Generated at Fri Jun 12 03:28:44 2026 by Isobel Mawby using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "larpandora/LArPandoraUtils/PandoraEventUtils.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/LArPandoraPIDData/PandoraPIDResult.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/TensorFlow/IvysaurusEvaluator.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/PFParticle.h"

#include <memory>

namespace pandorapid {
  class PandoraPID;
}

///////////////////////////////////////////////////////////////////////////////

class pandorapid::PandoraPID : public art::EDProducer {
public:
  explicit PandoraPID(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PandoraPID(PandoraPID const&) = delete;
  PandoraPID(PandoraPID&&) = delete;
  PandoraPID& operator=(PandoraPID const&) = delete;
  PandoraPID& operator=(PandoraPID&&) = delete;

  // Required functions.
  void produce(art::Event& evt) override;

private:

  // Module labels
  std::string m_recoModuleLabel;    

  // Algs 
  ivysaurus::IvysaurusEvaluator m_ivysaurusEvaluator;
};

///////////////////////////////////////////////////////////////////////////////

pandorapid::PandoraPID::PandoraPID(fhicl::ParameterSet const& pset) :
    EDProducer(pset),
    m_recoModuleLabel(pset.get<std::string>("RecoModuleLabel")),
    m_ivysaurusEvaluator(pset.get<fhicl::ParameterSet>("IvysaurusEvaluator"))      
{
    produces<std::vector<PandoraPIDResult>>();
    produces<art::Assns<recob::PFParticle, PandoraPIDResult, void>>();
}

///////////////////////////////////////////////////////////////////////////////

void pandorapid::PandoraPID::produce(art::Event& evt)
{
    std::unique_ptr<std::vector<PandoraPIDResult>> pandoraPIDResults(new std::vector<PandoraPIDResult>);
    std::unique_ptr<art::Assns<recob::PFParticle, PandoraPIDResult>> pfpAssns(new art::Assns<recob::PFParticle, PandoraPIDResult>);

    std::vector<art::Ptr<recob::PFParticle>> pfps = lar_pandora::PandoraEventUtils::GetPFParticles(evt, m_recoModuleLabel);
    for (const art::Ptr<recob::PFParticle> pfp : pfps)
    {
        // Ivysaurus
        ivysaurus::IvysaurusEvaluator::IvysaurusScores ivyScores = m_ivysaurusEvaluator.IvysaurusUseEvaluate(evt, pfp);

        PandoraPIDResult pandoraPIDResult;
        pandoraPIDResult.m_ivysaurusMuonScore = ivyScores.m_muonScore;
        pandoraPIDResult.m_ivysaurusProtonScore = ivyScores.m_protonScore;
        pandoraPIDResult.m_ivysaurusPionScore = ivyScores.m_pionScore;
        pandoraPIDResult.m_ivysaurusElectronScore = ivyScores.m_electronScore;
        pandoraPIDResult.m_ivysaurusPhotonScore = ivyScores.m_photonScore;
        pandoraPIDResult.m_ivysaurusParticleType = ivyScores.m_particleType;

        std::cout << "muon: " << ivyScores.m_muonScore;


        
        pandoraPIDResults->push_back(pandoraPIDResult);
        util::CreateAssn(evt, *pandoraPIDResults, pfp, *pfpAssns);
    }

    evt.put(std::move(pandoraPIDResults));
    evt.put(std::move(pfpAssns));
}

///////////////////////////////////////////////////////////////////////////////

DEFINE_ART_MODULE(pandorapid::PandoraPID)
