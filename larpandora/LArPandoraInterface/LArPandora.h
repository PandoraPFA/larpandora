/**
 *  @file   larpandora/LArPandoraInterface/LArPandora.h
 *
 *  @brief  Base producer module for reconstructing recob::PFParticles from recob::Hits
 *
 */

#ifndef LAR_PANDORA_H
#define LAR_PANDORA_H 1

#include "larpandora/LArPandoraInterface/ILArPandora.h"
#include "larpandora/LArPandoraInterface/LArPandoraInput.h"
#include "larpandora/LArPandoraInterface/LArPandoraOutput.h"

#include "larpandora/LArPandoraInterface/LArPandoraHitCollectionTool.h"

#include <string>

namespace lar_pandora {

  /**
 *  @brief  LArPandora class
 */
  class LArPandora : public ILArPandora {
  public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset the parameter set
     */
    LArPandora(fhicl::ParameterSet const& pset);

    void beginJob();
    void produce(art::Event& evt);

  protected:
    void CreatePandoraInput(art::Event& evt, IdToHitMap& idToHitMap);
    void ProcessPandoraOutput(art::Event& evt, const IdToHitMap& idToHitMap);

    fhicl::ParameterSet ConstructHitCollectionToolParameterSet(const fhicl::ParameterSet& pset);

    std::string m_configFile; ///< The config file

    bool
      m_shouldRunAllHitsCosmicReco; ///< Steering: whether to run all hits cosmic-ray reconstruction
    bool
      m_shouldRunStitching; ///< Steering: whether to stitch cosmic-ray muons crossing between volumes
    bool m_shouldRunCosmicHitRemoval; ///< Steering: whether to remove hits from tagged cosmic-rays
    bool
      m_shouldRunSlicing; ///< Steering: whether to slice events into separate regions for processing
    bool
      m_shouldRunNeutrinoRecoOption; ///< Steering: whether to run neutrino reconstruction for each slice
    bool
      m_shouldRunCosmicRecoOption; ///< Steering: whether to run cosmic-ray reconstruction for each slice
    bool
      m_shouldPerformSliceId; ///< Steering: whether to identify slices and select most appropriate pfos
    bool m_shouldProduceAllOutcomes; ///< Steering: whether to produce all reconstruction outcomes
    bool m_shouldCollectHitPredictions; ///< Whether there are some hit-level predicted labels to be collected
    bool m_printOverallRecoStatus; ///< Steering: whether to print current operation status messages

    std::string m_generatorModuleLabel;   ///< The generator module label
    std::string m_geantModuleLabel;       ///< The geant module label
    std::string m_simChannelModuleLabel;  ///< The SimChannel producer module label
    std::string m_hitfinderModuleLabel;   ///< The hit finder module label
    std::string m_backtrackerModuleLabel; ///< The back tracker module label

    std::string m_allOutcomesInstanceLabel; ///< The instance label for all outcomes

    bool m_enableProduction;   ///< Whether to persist output products
    bool m_enableDetectorGaps; ///< Whether to pass detector gap information to Pandora instances
    bool
      m_enableMCParticles; ///< Whether to pass mc information to Pandora instances to aid development
    bool
      m_disableRealDataCheck; ///< Whether to check if the input file contains real data before accessing MC information
    bool m_lineGapsCreated; ///< Book-keeping: whether line gap creation has been called

    std::unique_ptr<IHitCollectionTool> m_collectHitsTool; ///< art tool used to collect the hits

    LArPandoraInput::Settings m_inputSettings;   ///< The lar pandora input settings
    LArPandoraOutput::Settings m_outputSettings; ///< The lar pandora output settings

    LArDriftVolumeMap m_driftVolumeMap; ///< The map from volume id to drift volume
  };

} // namespace lar_pandora

#endif // #ifndef LAR_PANDORA_H
