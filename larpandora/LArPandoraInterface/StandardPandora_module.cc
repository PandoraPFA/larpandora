/**
 *  @file   larpandora/LArPandoraInterface/StandardPandora_module.cc
 *
 *  @brief  A generic LArPandora ART Producer module intended to work on ALL LAr-TPC wire-readout experiments.
 */

#include "art/Framework/Core/ModuleMacros.h"

#include "larpandora/LArPandoraInterface/LArPandora.h"

#include <string>

namespace lar_pandora
{

/**
 *  @brief  StandardPandora class
 */
class StandardPandora : public LArPandora
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset the parameter set
     */
    StandardPandora(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~StandardPandora();

private:
    void CreatePandoraInstances();
    void ConfigurePandoraInstances();
    void RunPandoraInstances();
    void ResetPandoraInstances();
    void DeletePandoraInstances();

    /**
     *  @brief  Pass external steering parameters, read from fhicl parameter set, to LArMaster Pandora algorithm
     *
     *  @param  pPandora the address of the relevant pandora instance
     */
    void ProvideExternalSteeringParameters(const pandora::Pandora *const pPandora) const;
};

DEFINE_ART_MODULE(StandardPandora)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "cetlib_except/exception.h"

#include "Api/PandoraApi.h"

#include "larpandoracontent/LArContent.h"
#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"
#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

#ifdef LIBTORCH_DL
    #include "larpandoradlcontent/LArDLContent.h"
#endif

namespace lar_pandora
{

StandardPandora::StandardPandora(fhicl::ParameterSet const &pset) :
    LArPandora(pset)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StandardPandora::~StandardPandora()
{
    this->DeletePandoraInstances();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StandardPandora::CreatePandoraInstances()
{
    m_pPrimaryPandora = new pandora::Pandora();
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*m_pPrimaryPandora));
#ifdef LIBTORCH_DL
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArDLContent::RegisterAlgorithms(*m_pPrimaryPandora));
#endif
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::RegisterBasicPlugins(*m_pPrimaryPandora));

    // ATTN Potentially ill defined, unless coordinate system set up to ensure that all drift volumes have same wire angles and pitches
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*m_pPrimaryPandora, new lar_content::LArPseudoLayerPlugin));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetLArTransformationPlugin(*m_pPrimaryPandora, new lar_content::LArRotationalTransformationPlugin));

    MultiPandoraApi::AddPrimaryPandoraInstance(m_pPrimaryPandora);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StandardPandora::ConfigurePandoraInstances()
{
    cet::search_path sp("FW_SEARCH_PATH");
    std::string fullConfigFileName;

    if (!sp.find_file(m_configFile, fullConfigFileName))
        throw cet::exception("StandardPandora") << " ConfigurePrimaryPandoraInstance - Failed to find xml configuration file " << m_configFile << " in FW search path";

    this->ProvideExternalSteeringParameters(m_pPrimaryPandora);
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*m_pPrimaryPandora, fullConfigFileName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StandardPandora::RunPandoraInstances()
{
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pPrimaryPandora));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StandardPandora::ResetPandoraInstances()
{
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pPrimaryPandora));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StandardPandora::DeletePandoraInstances()
{
    MultiPandoraApi::DeletePandoraInstances(m_pPrimaryPandora);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StandardPandora::ProvideExternalSteeringParameters(const pandora::Pandora *const pPandora) const
{
    auto *const pEventSteeringParameters = new lar_content::MasterAlgorithm::ExternalSteeringParameters;
    pEventSteeringParameters->m_shouldRunAllHitsCosmicReco = m_shouldRunAllHitsCosmicReco;
    pEventSteeringParameters->m_shouldRunStitching = m_shouldRunStitching;
    pEventSteeringParameters->m_shouldRunCosmicHitRemoval = m_shouldRunCosmicHitRemoval;
    pEventSteeringParameters->m_shouldRunSlicing = m_shouldRunSlicing;
    pEventSteeringParameters->m_shouldRunNeutrinoRecoOption = m_shouldRunNeutrinoRecoOption;
    pEventSteeringParameters->m_shouldRunCosmicRecoOption = m_shouldRunCosmicRecoOption;
    pEventSteeringParameters->m_shouldPerformSliceId = m_shouldPerformSliceId;
    pEventSteeringParameters->m_printOverallRecoStatus = m_printOverallRecoStatus;
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, pandora::ExternallyConfiguredAlgorithm::SetExternalParameters(*pPandora, "LArMaster", pEventSteeringParameters));

#ifdef LIBTORCH_DL
    auto *const pEventSettingsParametersCopy = new lar_content::MasterAlgorithm::ExternalSteeringParameters(*pEventSteeringParameters);
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, pandora::ExternallyConfiguredAlgorithm::SetExternalParameters(*pPandora,
        "LArDLMaster", pEventSettingsParametersCopy));
#endif
}

} // namespace lar_pandora
