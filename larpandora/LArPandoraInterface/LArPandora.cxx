/**
 *  @file   larpandora/LArPandoraInterface/LArPandora.cxx
 *
 *  @brief  Base producer module for reconstructing recob::PFParticles from recob::Hits
 *
 */

#include "larpandora/LArPandoraInterface/LArPandora.h"

#include "art/Framework/Principal/Event.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "cetlib/cpu_timer.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "larcore/Geometry/Geometry.h"

#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include "Pandora/PandoraInputTypes.h"
#include "Api/PandoraApi.h"

#include "larpandoracontent/LArContent.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include <iostream>
#include <limits>

namespace lar_pandora {

  LArPandora::LArPandora(fhicl::ParameterSet const& pset)
    : ILArPandora(pset)
    , m_configFile(pset.get<std::string>("ConfigFile"))
    , m_shouldRunAllHitsCosmicReco(pset.get<bool>("ShouldRunAllHitsCosmicReco"))
    , m_shouldRunStitching(pset.get<bool>("ShouldRunStitching"))
    , m_shouldRunCosmicHitRemoval(pset.get<bool>("ShouldRunCosmicHitRemoval"))
    , m_shouldRunSlicing(pset.get<bool>("ShouldRunSlicing"))
    , m_shouldRunNeutrinoRecoOption(pset.get<bool>("ShouldRunNeutrinoRecoOption"))
    , m_shouldRunCosmicRecoOption(pset.get<bool>("ShouldRunCosmicRecoOption"))
    , m_shouldPerformSliceId(pset.get<bool>("ShouldPerformSliceId"))
    , m_shouldProduceAllOutcomes(pset.get<bool>("ProduceAllOutcomes", false))
    , m_printOverallRecoStatus(pset.get<bool>("PrintOverallRecoStatus", false))
    , m_generatorModuleLabel(pset.get<std::string>("GeneratorModuleLabel", ""))
    , m_geantModuleLabel(pset.get<std::string>("GeantModuleLabel", "largeant"))
    , m_simChannelModuleLabel(pset.get<std::string>("SimChannelModuleLabel", m_geantModuleLabel))
    , m_hitfinderModuleLabel(pset.get<std::string>("HitFinderModuleLabel"))
    , m_backtrackerModuleLabel(pset.get<std::string>("BackTrackerModuleLabel", ""))
    , m_allOutcomesInstanceLabel(pset.get<std::string>("AllOutcomesInstanceLabel", "allOutcomes"))
    , m_enableProduction(pset.get<bool>("EnableProduction", true))
    , m_enableDetectorGaps(pset.get<bool>("EnableLineGaps", true))
    , m_enableMCParticles(pset.get<bool>("EnableMCParticles", false))
    , m_disableRealDataCheck(pset.get<bool>("DisableRealDataCheck", false))
    , m_lineGapsCreated(false)
    , m_collectHitsTool{
        art::make_tool<IHitCollectionTool>(this->ConstructHitCollectionToolParameterSet(pset))}
  {
    m_inputSettings.m_useHitWidths = pset.get<bool>("UseHitWidths", true);
    m_inputSettings.m_useBirksCorrection = pset.get<bool>("UseBirksCorrection", false);
    m_inputSettings.m_useActiveBoundingBox = pset.get<bool>("UseActiveBoundingBox", false);
    m_inputSettings.m_uidOffset = pset.get<int>("UidOffset", 100000000);
    m_inputSettings.m_dx_cm = pset.get<double>("DefaultHitWidth", 0.5);
    m_inputSettings.m_int_cm = pset.get<double>("InteractionLength", 84.);
    m_inputSettings.m_rad_cm = pset.get<double>("RadiationLength", 14.);
    m_inputSettings.m_dEdX_mip = pset.get<double>("dEdXmip", 2.);
    m_inputSettings.m_mips_max = pset.get<double>("MipsMax", 50.);
    m_inputSettings.m_mips_if_negative = pset.get<double>("MipsIfNegative", 0.);
    m_inputSettings.m_mips_to_gev = pset.get<double>("MipsToGeV", 3.5e-4);
    m_inputSettings.m_recombination_factor = pset.get<double>("RecombinationFactor", 0.63);
    m_outputSettings.m_shouldRunStitching = m_shouldRunStitching;
    m_outputSettings.m_shouldProduceSlices = pset.get<bool>("ShouldProduceSlices", true);
    m_outputSettings.m_shouldProduceTestBeamInteractionVertices =
      pset.get<bool>("ShouldProduceTestBeamInteractionVertices", false);
    m_outputSettings.m_testBeamInteractionVerticesInstanceLabel = pset.get<std::string>(
      "TestBeamInteractionVerticesInstanceLabel", "testBeamInteractionVertices");
    m_outputSettings.m_isNeutrinoRecoOnlyNoSlicing =
      (!m_shouldRunSlicing && m_shouldRunNeutrinoRecoOption && !m_shouldRunCosmicRecoOption);
    m_outputSettings.m_hitfinderModuleLabel = m_hitfinderModuleLabel;

    if (m_enableProduction) {
      // Set up the instance names to produces
      std::vector<std::string> instanceNames({""});
      if (m_shouldProduceAllOutcomes) instanceNames.push_back(m_allOutcomesInstanceLabel);

      for (const std::string& instanceName : instanceNames) {
        produces<std::vector<recob::PFParticle>>(instanceName);
        produces<std::vector<recob::SpacePoint>>(instanceName);
        produces<std::vector<recob::Cluster>>(instanceName);
        produces<std::vector<recob::Vertex>>(instanceName);
        produces<std::vector<larpandoraobj::PFParticleMetadata>>(instanceName);

        produces<art::Assns<recob::PFParticle, larpandoraobj::PFParticleMetadata>>(instanceName);
        produces<art::Assns<recob::PFParticle, recob::SpacePoint>>(instanceName);
        produces<art::Assns<recob::PFParticle, recob::Cluster>>(instanceName);
        produces<art::Assns<recob::PFParticle, recob::Vertex>>(instanceName);
        produces<art::Assns<recob::SpacePoint, recob::Hit>>(instanceName);
        produces<art::Assns<recob::Cluster, recob::Hit>>(instanceName);

        if (m_outputSettings.m_shouldProduceTestBeamInteractionVertices) {
          // ATTN: Test beam interaction vertex instance label appended to current instance name to preserve unique label in multiple instance case
          produces<std::vector<recob::Vertex>>(
            instanceName + m_outputSettings.m_testBeamInteractionVerticesInstanceLabel);
          produces<art::Assns<recob::PFParticle, recob::Vertex>>(
            instanceName + m_outputSettings.m_testBeamInteractionVerticesInstanceLabel);
        }

        if (m_outputSettings.m_shouldRunStitching) {
          produces<std::vector<anab::T0>>(instanceName);
          produces<art::Assns<recob::PFParticle, anab::T0>>(instanceName);
        }

        if (m_outputSettings.m_shouldProduceSlices) {
          produces<std::vector<recob::Slice>>(instanceName);
          produces<art::Assns<recob::Slice, recob::Hit>>(instanceName);
          produces<art::Assns<recob::PFParticle, recob::Slice>>(instanceName);
        }
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandora::beginJob()
  {
    LArDriftVolumeList driftVolumeList;
    LArPandoraGeometry::LoadGeometry(
      driftVolumeList, m_driftVolumeMap, m_inputSettings.m_useActiveBoundingBox);

    this->CreatePandoraInstances();

    if (!m_pPrimaryPandora)
      throw cet::exception("LArPandora")
        << " LArPandora::beginJob - failed to create primary Pandora instance " << std::endl;

    m_inputSettings.m_pPrimaryPandora = m_pPrimaryPandora;
    m_outputSettings.m_pPrimaryPandora = m_pPrimaryPandora;

    // Pass basic LArTPC information to pandora instances
    LArPandoraInput::CreatePandoraLArTPCs(m_inputSettings, driftVolumeList);

    // If using global drift volume approach, pass details of gaps between daughter volumes to the pandora instance
    if (m_enableDetectorGaps) {
      LArDetectorGapList listOfGaps;
      LArPandoraGeometry::LoadDetectorGaps(listOfGaps, m_inputSettings.m_useActiveBoundingBox);
      LArPandoraInput::CreatePandoraDetectorGaps(m_inputSettings, driftVolumeList, listOfGaps);
    }

    // Parse Pandora settings xml files
    this->ConfigurePandoraInstances();
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandora::produce(art::Event& evt)
  {
    this->SetPandoraEventInformation(evt);
    IdToHitMap idToHitMap;
    this->CreatePandoraInput(evt, idToHitMap);
    this->RunPandoraInstances();
    this->ProcessPandoraOutput(evt, idToHitMap);
    this->ResetPandoraInstances();
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandora::SetPandoraEventInformation(art::Event& evt)
  {

    try
    {
        pandora::PandoraInputType run(evt.run());
        pandora::PandoraInputType subrun(evt.subRun());
        pandora::PandoraInputType event(evt.event());
        PandoraApi::SetEventInformation(*m_pPrimaryPandora, run, subrun, event);
    }
    catch (...)
    {
        throw cet::exception("LArPandora") << "LArPandora::SetPandoraEventInformation --- unable to "
            "form PandoraInputType from R:S:E info";
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------


  void LArPandora::CreatePandoraInput(art::Event& evt, IdToHitMap& idToHitMap)
  {
    // ATTN Should complete gap creation in begin job callback, but channel status service functionality unavailable at that point
    if (!m_lineGapsCreated && m_enableDetectorGaps) {
      LArPandoraInput::CreatePandoraReadoutGaps(m_inputSettings, m_driftVolumeMap);
      m_lineGapsCreated = true;
    }

    HitVector artHits;
    SimChannelVector artSimChannels;
    HitsToTrackIDEs artHitsToTrackIDEs;
    MCParticleVector artMCParticleVector;
    RawMCParticleVector generatorArtMCParticleVector;
    MCTruthToMCParticles artMCTruthToMCParticles;
    MCParticlesToMCTruth artMCParticlesToMCTruth;

    bool areSimChannelsValid(false);

    m_collectHitsTool->CollectHits(evt, m_hitfinderModuleLabel, artHits);

    if (m_enableMCParticles && (m_disableRealDataCheck || !evt.isRealData())) {
      LArPandoraHelper::CollectMCParticles(evt, m_geantModuleLabel, artMCParticleVector);

      if (!m_generatorModuleLabel.empty())
        LArPandoraHelper::CollectGeneratorMCParticles(
          evt, m_generatorModuleLabel, generatorArtMCParticleVector);

      LArPandoraHelper::CollectMCParticles(
        evt, m_geantModuleLabel, artMCTruthToMCParticles, artMCParticlesToMCTruth);

      LArPandoraHelper::CollectSimChannels(
        evt, m_simChannelModuleLabel, artSimChannels, areSimChannelsValid);
      if (!artSimChannels.empty()) {
        LArPandoraHelper::BuildMCParticleHitMaps(evt, artHits, artSimChannels, artHitsToTrackIDEs);
      }
      else if (!areSimChannelsValid) {
        if (m_backtrackerModuleLabel.empty())
          throw cet::exception("LArPandora")
            << "LArPandora::CreatePandoraInput - Can't build MCParticle to Hit map." << std::endl
            << "No SimChannels found with label \"" << m_simChannelModuleLabel
            << "\", and BackTrackerModuleLabel isn't set in FHiCL." << std::endl;

        LArPandoraHelper::BuildMCParticleHitMaps(
          evt, m_hitfinderModuleLabel, m_backtrackerModuleLabel, artHitsToTrackIDEs);
      }
      else {
        mf::LogDebug("LArPandora")
          << " *** LArPandora::CreatePandoraInput - empty list of sim channels found " << std::endl;
      }
    }

    LArPandoraInput::CreatePandoraHits2D(
      evt, m_inputSettings, m_driftVolumeMap, artHits, idToHitMap);

    if (m_enableMCParticles && (m_disableRealDataCheck || !evt.isRealData())) {
      LArPandoraInput::CreatePandoraMCParticles(m_inputSettings,
                                                artMCTruthToMCParticles,
                                                artMCParticlesToMCTruth,
                                                generatorArtMCParticleVector);
      LArPandoraInput::CreatePandoraMCLinks2D(m_inputSettings, idToHitMap, artHitsToTrackIDEs);
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandora::ProcessPandoraOutput(art::Event& evt, const IdToHitMap& idToHitMap)
  {
    if (m_enableProduction) {
      m_outputSettings.m_shouldProduceAllOutcomes = false;
      LArPandoraOutput::ProduceArtOutput(m_outputSettings, idToHitMap, evt);

      if (m_shouldProduceAllOutcomes) {
        m_outputSettings.m_shouldProduceAllOutcomes = true;
        m_outputSettings.m_allOutcomesInstanceLabel = m_allOutcomesInstanceLabel;
        LArPandoraOutput::ProduceArtOutput(m_outputSettings, idToHitMap, evt);
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  fhicl::ParameterSet LArPandora::ConstructHitCollectionToolParameterSet(
    const fhicl::ParameterSet& pset)
  {
    if (pset.has_key("HitCollectionTool")) {
      return pset.get<fhicl::ParameterSet>("HitCollectionTool");
    }

    fhicl::ParameterSet psetHitCollection;
    psetHitCollection.put<std::string>("tool_type", "LArPandoraHitCollectionToolDefault");
    return psetHitCollection;
  }

} // namespace lar_pandora
