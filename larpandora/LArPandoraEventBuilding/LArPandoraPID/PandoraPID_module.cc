/**
 *  @file   larpandora/LArPandoraEventBuilding/LArPandoraPID/PandoraPID_module.cc
 *
 *  @brief  PandoraPIDResult producer
 */
// ART
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
// LArSoft
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraUtils/PandoraEventUtils.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/LArPandoraPIDData/PandoraPIDResult.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Evaluator/IvysaurusEvaluator.h"
// ROOT
//#include <memory>

namespace pandorapid {

class PandoraPID : public art::EDProducer {
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
  std::string m_recoModuleLabel;    
  ivysaurus::IvysaurusEvaluator m_ivysaurusEvaluator;
};

//------------------------------------------------------------------------------------------------------------------------------------------

PandoraPID::PandoraPID(fhicl::ParameterSet const& pset) :
    EDProducer(pset),
    m_recoModuleLabel(pset.get<std::string>("RecoModuleLabel")),
    m_ivysaurusEvaluator(pset.get<fhicl::ParameterSet>("IvysaurusEvaluator"))      
{
    produces<std::vector<PandoraPIDResult>>();
    produces<art::Assns<recob::PFParticle, PandoraPIDResult, void>>();
}


//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPID::produce(art::Event& evt)
{
    std::unique_ptr<std::vector<PandoraPIDResult>> pandoraPIDResults(new std::vector<PandoraPIDResult>);
    std::unique_ptr<art::Assns<recob::PFParticle, PandoraPIDResult>> pfpAssns(new art::Assns<recob::PFParticle, PandoraPIDResult>);

    std::vector<art::Ptr<recob::PFParticle>> pfps = lar_pandora::PandoraEventUtils::GetPFParticles(evt, m_recoModuleLabel);
    for (const art::Ptr<recob::PFParticle> pfp : pfps)
    {
        ivysaurus::IvysaurusEvaluator::IvysaurusScores ivyScores = m_ivysaurusEvaluator.IvysaurusUseEvaluate(evt, pfp);

        PandoraPIDResult pandoraPIDResult;
        pandoraPIDResult.SetIvysaurusScores(ivyScores.m_muonScore, ivyScores.m_protonScore, ivyScores.m_pionScore,
            ivyScores.m_electronScore, ivyScores.m_photonScore);  

        pandoraPIDResults->push_back(pandoraPIDResult);
        util::CreateAssn(evt, *pandoraPIDResults, pfp, *pfpAssns);
    }

    evt.put(std::move(pandoraPIDResults));
    evt.put(std::move(pfpAssns));
}

} // namespace pandorapid

DEFINE_ART_MODULE(pandorapid::PandoraPID)

