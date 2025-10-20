/**
 *  @file   larpandora/LArPandoraInterface/LArPandoraInput.cxx
 *
 *  @brief  Helper functions for providing inputs to pandora
 */

#include "larpandora/LArPandoraInterface/LArPandoraInput.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

#include "lardataobj/RecoBase/Hit.h"

#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"

#include "nusimdata/SimulationBase/MCTruth.h"

#include "larsim/MCCheater/ParticleInventoryService.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"

#include "Api/PandoraApi.h"
#include "Managers/PluginManager.h"
#include "Plugins/LArTransformationPlugin.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "larpandora/LArPandoraInterface/Detectors/GetDetectorType.h"
#include "larpandora/LArPandoraInterface/Detectors/LArPandoraDetectorType.h"
#include "larpandora/LArPandoraInterface/ILArPandora.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include <limits>
#include <utility>

namespace lar_pandora {

  void LArPandoraInput::CreatePandoraHits2D(const art::Event& e,
                                            const Settings& settings,
                                            const LArDriftVolumeMap& driftVolumeMap,
                                            const HitVector& hitVector,
                                            const HitToScores& hitToScores,
                                            const HitToScoreLabels& hitToScoreLabels,
                                            IdToHitMap& idToHitMap)
  {
    mf::LogDebug("LArPandora") << " *** LArPandoraInput::CreatePandoraHits2D(...) *** "
                               << std::endl;

    if (!settings.m_pPrimaryPandora)
      throw cet::exception("LArPandora")
        << "CreatePandoraHits2D - primary Pandora instance does not exist ";

    const pandora::Pandora* pPandora(settings.m_pPrimaryPandora);

    auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout const>()->Get();
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);
    LArPandoraDetectorType* detType(detector_functions::GetDetectorType());

    // Loop over ART hits
    int hitCounter(settings.m_hitCounterOffset);

    lar_content::LArCaloHitFactory caloHitFactory;

    for (auto const& hit : hitVector) {
      const geo::WireID hit_WireID(hit->WireID());

      // Get basic hit properties (view, time, charge)
      const geo::View_t hit_View(hit->View());
      const double hit_Charge(hit->Integral());
      const double hit_Time(hit->PeakTime());
      const double hit_TimeStart(hit->PeakTimeMinusRMS());
      const double hit_TimeEnd(hit->PeakTimePlusRMS());

      // Get hit X coordinate and, if using a single global drift volume, remove any out-of-time hits here
      const double xpos_cm(detProp.ConvertTicksToX(hit_Time, hit_WireID.parentID()));
      const double dxpos_cm(
        std::fabs(detProp.ConvertTicksToX(hit_TimeEnd, hit_WireID.parentID()) -
                  detProp.ConvertTicksToX(hit_TimeStart, hit_WireID.parentID())));

      // Get hit Y and Z coordinates, based on central position of wire
      auto const xyz = wireReadoutGeom.Wire(hit_WireID).GetCenter();
      const double y0_cm(xyz.Y());
      const double z0_cm(xyz.Z());

      // Get other hit properties here
      const double wire_pitch_cm(wireReadoutGeom.Plane({0, 0}, hit_View).WirePitch()); // cm
      const double mips(LArPandoraInput::GetMips(detProp, settings, hit_Charge, hit_View));

      // Create Pandora CaloHit
      lar_content::LArCaloHitParameters caloHitParameters;

      try {
        caloHitParameters.m_expectedDirection = pandora::CartesianVector(0., 0., 1.);
        caloHitParameters.m_cellNormalVector = pandora::CartesianVector(0., 0., 1.);
        caloHitParameters.m_cellSize0 = settings.m_dx_cm;
        caloHitParameters.m_cellSize1 = (settings.m_useHitWidths ? dxpos_cm : settings.m_dx_cm);
        caloHitParameters.m_cellThickness = wire_pitch_cm;
        caloHitParameters.m_cellGeometry = pandora::RECTANGULAR;
        caloHitParameters.m_time = 0.;
        caloHitParameters.m_nCellRadiationLengths = settings.m_dx_cm / settings.m_rad_cm;
        caloHitParameters.m_nCellInteractionLengths = settings.m_dx_cm / settings.m_int_cm;
        caloHitParameters.m_isDigital = false;
        caloHitParameters.m_hitRegion = pandora::SINGLE_REGION;
        caloHitParameters.m_layer = 0;
        caloHitParameters.m_isInOuterSamplingLayer = false;
        caloHitParameters.m_inputEnergy = hit_Charge;
        caloHitParameters.m_mipEquivalentEnergy = mips;
        caloHitParameters.m_electromagneticEnergy = mips * settings.m_mips_to_gev;
        caloHitParameters.m_hadronicEnergy = mips * settings.m_mips_to_gev;
        caloHitParameters.m_pParentAddress = (void*)((intptr_t)(++hitCounter));
        caloHitParameters.m_larTPCVolumeId =
          LArPandoraGeometry::GetVolumeID(driftVolumeMap, hit_WireID.Cryostat, hit_WireID.TPC);
        caloHitParameters.m_daughterVolumeId = LArPandoraGeometry::GetDaughterVolumeID(
          driftVolumeMap, hit_WireID.Cryostat, hit_WireID.TPC);

        if (hit_View == detType->TargetViewW(hit_WireID.TPC, hit_WireID.Cryostat)) {
          caloHitParameters.m_hitType = pandora::TPC_VIEW_W;
          const double wpos_cm(
            pPandora->GetPlugins()->GetLArTransformationPlugin()->YZtoW(y0_cm, z0_cm));
          caloHitParameters.m_positionVector = pandora::CartesianVector(xpos_cm, 0., wpos_cm);
        }
        else if (hit_View == detType->TargetViewU(hit_WireID.TPC, hit_WireID.Cryostat)) {
          caloHitParameters.m_hitType = pandora::TPC_VIEW_U;
          const double upos_cm(
            pPandora->GetPlugins()->GetLArTransformationPlugin()->YZtoU(y0_cm, z0_cm));
          caloHitParameters.m_positionVector = pandora::CartesianVector(xpos_cm, 0., upos_cm);
        }
        else if (hit_View == detType->TargetViewV(hit_WireID.TPC, hit_WireID.Cryostat)) {
          caloHitParameters.m_hitType = pandora::TPC_VIEW_V;
          const double vpos_cm(
            pPandora->GetPlugins()->GetLArTransformationPlugin()->YZtoV(y0_cm, z0_cm));
          caloHitParameters.m_positionVector = pandora::CartesianVector(xpos_cm, 0., vpos_cm);
        }
        else {
          throw cet::exception("LArPandora")
            << "CreatePandoraHits2D - this wire view not recognised (View=" << hit_View << ") ";
        }
      }
      catch (const pandora::StatusCodeException&) {
        mf::LogWarning("LArPandora")
          << "CreatePandoraHits2D - invalid calo hit parameter provided, all assigned values must "
             "be finite, calo hit omitted "
          << std::endl;
        continue;
      }

      // If wanted, populate the hit predictions
      if (settings.m_useHitPredictions) {
        auto itScorePred = hitToScores.find(hit);
        auto itLabelPred = hitToScoreLabels.find(hit);
        if ((itScorePred != hitToScores.end()) && 
            (itLabelPred != hitToScoreLabels.end())) {
          caloHitParameters.m_hitScores = itScorePred->second;
          caloHitParameters.m_hitScoreLabels = itLabelPred->second;
        } 
      }

      // Store the hit address
      if (hitCounter >= settings.m_uidOffset)
        throw cet::exception("LArPandora")
          << "CreatePandoraHits2D - detected an excessive number of hits (" << hitCounter << ") ";

      idToHitMap[hitCounter] = hit;

      // Create the Pandora hit
      try {
        PANDORA_THROW_RESULT_IF(
          pandora::STATUS_CODE_SUCCESS,
          !=,
          PandoraApi::CaloHit::Create(*pPandora, caloHitParameters, caloHitFactory));
      }
      catch (const pandora::StatusCodeException&) {
        mf::LogWarning("LArPandora") << "CreatePandoraHits2D - unable to create calo hit, "
                                        "insufficient or invalid information supplied "
                                     << std::endl;
        continue;
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraInput::CreatePandoraLArTPCs(const Settings& settings,
                                             const LArDriftVolumeList& driftVolumeList)
  {
    mf::LogDebug("LArPandora") << " *** LArPandoraInput::CreatePandoraLArTPCs(...) *** "
                               << std::endl;

    if (!settings.m_pPrimaryPandora)
      throw cet::exception("LArPandora")
        << "CreatePandoraDetectorGaps - primary Pandora instance does not exist ";

    const pandora::Pandora* pPandora(settings.m_pPrimaryPandora);

    for (const LArDriftVolume& driftVolume : driftVolumeList) {
      PandoraApi::Geometry::LArTPC::Parameters parameters;

      try {
        parameters.m_larTPCVolumeId = driftVolume.GetVolumeID();
        parameters.m_centerX = driftVolume.GetCenterX();
        parameters.m_centerY = driftVolume.GetCenterY();
        parameters.m_centerZ = driftVolume.GetCenterZ();
        parameters.m_widthX = driftVolume.GetWidthX();
        parameters.m_widthY = driftVolume.GetWidthY();
        parameters.m_widthZ = driftVolume.GetWidthZ();
        parameters.m_wirePitchU = driftVolume.GetWirePitchU();
        parameters.m_wirePitchV = driftVolume.GetWirePitchV();
        parameters.m_wirePitchW = driftVolume.GetWirePitchW();
        parameters.m_wireAngleU = driftVolume.GetWireAngleU();
        parameters.m_wireAngleV = driftVolume.GetWireAngleV();
        parameters.m_wireAngleW = driftVolume.GetWireAngleW();
        parameters.m_sigmaUVW = driftVolume.GetSigmaUVZ();
        parameters.m_isDriftInPositiveX = driftVolume.IsPositiveDrift();
      }
      catch (const pandora::StatusCodeException&) {
        mf::LogWarning("LArPandora") << "CreatePandoraLArTPCs - invalid tpc parameter provided, "
                                        "all assigned values must be finite, tpc omitted "
                                     << std::endl;
        continue;
      }

      try {
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS,
                                !=,
                                PandoraApi::Geometry::LArTPC::Create(*pPandora, parameters));
      }
      catch (const pandora::StatusCodeException&) {
        mf::LogWarning("LArPandora") << "CreatePandoraLArTPCs - unable to create tpc, insufficient "
                                        "or invalid information supplied "
                                     << std::endl;
        continue;
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraInput::CreatePandoraDetectorGaps(const Settings& settings,
                                                  const LArDriftVolumeList& driftVolumeList,
                                                  const LArDetectorGapList& listOfGaps)
  {
    //ATTN - Unlike SP, DP detector gaps are not in the drift direction
    art::ServiceHandle<geo::Geometry const> theGeometry;
    LArPandoraDetectorType* detType(detector_functions::GetDetectorType());

    mf::LogDebug("LArPandora") << " *** LArPandoraInput::CreatePandoraDetectorGaps(...) *** "
                               << std::endl;

    if (!settings.m_pPrimaryPandora)
      throw cet::exception("LArPandora")
        << "CreatePandoraDetectorGaps - primary Pandora instance does not exist ";

    const pandora::Pandora* pPandora(settings.m_pPrimaryPandora);

    for (const LArDetectorGap& gap : listOfGaps) {
      PandoraApi::Geometry::LineGap::Parameters parameters;

      try {
        parameters = detType->CreateLineGapParametersFromDetectorGaps(gap);
      }
      catch (const pandora::StatusCodeException&) {
        mf::LogWarning("LArPandora")
          << "CreatePandoraDetectorGaps - invalid line gap parameter provided, all assigned values "
             "must be finite, line gap omitted "
          << std::endl;
        continue;
      }
      try {
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS,
                                !=,
                                PandoraApi::Geometry::LineGap::Create(*pPandora, parameters));
      }
      catch (const pandora::StatusCodeException&) {
        mf::LogWarning("LArPandora") << "CreatePandoraDetectorGaps - unable to create line gap, "
                                        "insufficient or invalid information supplied "
                                     << std::endl;
        continue;
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraInput::CreatePandoraReadoutGaps(const Settings& settings,
                                                 const LArDriftVolumeMap& driftVolumeMap)
  {
    mf::LogDebug("LArPandora") << " *** LArPandoraInput::CreatePandoraReadoutGaps(...) *** "
                               << std::endl;

    if (!settings.m_pPrimaryPandora)
      throw cet::exception("LArPandora")
        << "CreatePandoraReadoutGaps - primary Pandora instance does not exist ";

    const pandora::Pandora* pPandora(settings.m_pPrimaryPandora);

    const lariov::ChannelStatusProvider& channelStatus(
      art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider());

    LArPandoraDetectorType* detType(detector_functions::GetDetectorType());

    auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout const>()->Get();
    for (auto const& plane : wireReadoutGeom.Iterate<geo::PlaneGeo>()) {
      const float halfWirePitch(0.5f * plane.WirePitch());
      const unsigned int nWires(wireReadoutGeom.Nwires(plane.ID()));

      int firstBadWire(-1), lastBadWire(-1);

      for (unsigned int iwire = 0; iwire < nWires; ++iwire) {
        const raw::ChannelID_t channel(
          wireReadoutGeom.PlaneWireToChannel(geo::WireID{plane.ID(), iwire}));
        const bool isBadChannel(channelStatus.IsBad(channel));
        const bool isLastWire(nWires == (iwire + 1));

        if (isBadChannel && (firstBadWire < 0)) firstBadWire = iwire;

        if (isBadChannel || isLastWire) lastBadWire = iwire;

        if (isBadChannel && !isLastWire) continue;

        if ((firstBadWire < 0) || (lastBadWire < 0)) continue;

        auto const firstXYZ = plane.Wire(firstBadWire).GetCenter();
        auto const lastXYZ = plane.Wire(lastBadWire).GetCenter();

        firstBadWire = -1;
        lastBadWire = -1;

        PandoraApi::Geometry::LineGap::Parameters parameters;

        try {
          float xFirst(-std::numeric_limits<float>::max());
          float xLast(std::numeric_limits<float>::max());

          auto const [icstat, itpc] = std::make_pair(plane.ID().Cryostat, plane.ID().TPC);
          const unsigned int volumeId(
            LArPandoraGeometry::GetVolumeID(driftVolumeMap, icstat, itpc));
          LArDriftVolumeMap::const_iterator volumeIter(driftVolumeMap.find(volumeId));

          if (driftVolumeMap.end() != volumeIter) {
            xFirst = volumeIter->second.GetCenterX() - 0.5f * volumeIter->second.GetWidthX();
            xLast = volumeIter->second.GetCenterX() + 0.5f * volumeIter->second.GetWidthX();
          }

          const geo::View_t iview = plane.View();
          parameters = detType->CreateLineGapParametersFromReadoutGaps(
            iview, itpc, icstat, firstXYZ, lastXYZ, halfWirePitch, xFirst, xLast, pPandora);
        }
        catch (const pandora::StatusCodeException&) {
          mf::LogWarning("LArPandora")
            << "CreatePandoraReadoutGaps - invalid line gap parameter provided, all assigned "
               "values must be finite, line gap omitted "
            << std::endl;
          continue;
        }

        try {
          PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS,
                                  !=,
                                  PandoraApi::Geometry::LineGap::Create(*pPandora, parameters));
        }
        catch (const pandora::StatusCodeException&) {
          mf::LogWarning("LArPandora") << "CreatePandoraReadoutGaps - unable to create line "
                                          "gap, insufficient or invalid information supplied "
                                       << std::endl;
          continue;
        }
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraInput::CreatePandoraMCParticles(
    const Settings& settings,
    const MCTruthToMCParticles& truthToParticleMap,
    const MCParticlesToMCTruth& particleToTruthMap,
    const RawMCParticleVector& generatorMCParticleVector)
  {
    mf::LogDebug("LArPandora") << " *** LArPandoraInput::CreatePandoraMCParticles(...) *** "
                               << std::endl;
    art::ServiceHandle<cheat::ParticleInventoryService const> particleInventoryService;

    if (!settings.m_pPrimaryPandora)
      throw cet::exception("LArPandora")
        << "CreatePandoraMCParticles - primary Pandora instance does not exist ";

    const pandora::Pandora* pPandora(settings.m_pPrimaryPandora);

    // Make indexed list of MC particles
    MCParticleMap particleMap;

    for (MCParticlesToMCTruth::const_iterator iter = particleToTruthMap.begin(),
                                              iterEnd = particleToTruthMap.end();
         iter != iterEnd;
         ++iter) {
      const art::Ptr<simb::MCParticle> particle = iter->first;
      particleMap[particle->TrackId()] = particle;
    }

    // Loop over MC truth objects
    int neutrinoCounter(0);

    lar_content::LArMCParticleFactory mcParticleFactory;

    for (MCTruthToMCParticles::const_iterator iter1 = truthToParticleMap.begin(),
                                              iterEnd1 = truthToParticleMap.end();
         iter1 != iterEnd1;
         ++iter1) {
      const art::Ptr<simb::MCTruth> truth = iter1->first;

      if (truth->NeutrinoSet()) {
        const simb::MCNeutrino neutrino(truth->GetNeutrino());
        ++neutrinoCounter;

        if (neutrinoCounter >= settings.m_uidOffset)
          throw cet::exception("LArPandora")
            << "CreatePandoraMCParticles - detected an excessive number of mc neutrinos ("
            << neutrinoCounter << ")";

        const int neutrinoID(neutrinoCounter + 9 * settings.m_uidOffset);

        // Create Pandora 3D MC Particle
        lar_content::LArMCParticleParameters mcParticleParameters;

        try {
          if (truth->Origin() == simb::kSuperNovaNeutrino)
            mcParticleParameters.m_nuanceCode = 4000;
          else
            mcParticleParameters.m_nuanceCode = neutrino.InteractionType();
          mcParticleParameters.m_process = lar_content::MC_PROC_INCIDENT_NU;
          mcParticleParameters.m_energy = neutrino.Nu().E();
          mcParticleParameters.m_momentum =
            pandora::CartesianVector(neutrino.Nu().Px(), neutrino.Nu().Py(), neutrino.Nu().Pz());
          mcParticleParameters.m_vertex =
            pandora::CartesianVector(neutrino.Nu().Vx(), neutrino.Nu().Vy(), neutrino.Nu().Vz());
          mcParticleParameters.m_endpoint =
            pandora::CartesianVector(neutrino.Nu().Vx(), neutrino.Nu().Vy(), neutrino.Nu().Vz());
          mcParticleParameters.m_particleId = neutrino.Nu().PdgCode();
          mcParticleParameters.m_mcParticleType = pandora::MC_3D;
          mcParticleParameters.m_pParentAddress = (void*)((intptr_t)neutrinoID);
        }
        catch (const pandora::StatusCodeException&) {
          mf::LogWarning("LArPandora")
            << "CreatePandoraMCParticles - invalid mc neutrino parameter provided, all assigned "
               "values must be finite, mc neutrino omitted "
            << std::endl;
          continue;
        }

        try {
          PANDORA_THROW_RESULT_IF(
            pandora::STATUS_CODE_SUCCESS,
            !=,
            PandoraApi::MCParticle::Create(*pPandora, mcParticleParameters, mcParticleFactory));
        }
        catch (const pandora::StatusCodeException&) {
          mf::LogWarning("LArPandora") << "CreatePandoraMCParticles - unable to create mc "
                                          "neutrino, insufficient or invalid information supplied "
                                       << std::endl;
          continue;
        }

        // Loop over associated particles
        const MCParticleVector& particleVector = iter1->second;

        for (MCParticleVector::const_iterator iter2 = particleVector.begin(),
                                              iterEnd2 = particleVector.end();
             iter2 != iterEnd2;
             ++iter2) {
          const art::Ptr<simb::MCParticle> particle = *iter2;
          const int trackID(particle->TrackId());

          // Mother/Daughter Links
          if (particle->Mother() == 0) {
            try {
              PANDORA_THROW_RESULT_IF(
                pandora::STATUS_CODE_SUCCESS,
                !=,
                PandoraApi::SetMCParentDaughterRelationship(
                  *pPandora, (void*)((intptr_t)neutrinoID), (void*)((intptr_t)trackID)));
            }
            catch (const pandora::StatusCodeException&) {
              mf::LogWarning("LArPandora") << "CreatePandoraMCParticles - unable to create mc "
                                              "particle relationship, invalid information supplied "
                                           << std::endl;
              continue;
            }
          }
        }
      }
    }

    mf::LogDebug("LArPandora") << "   Number of Pandora neutrinos: " << neutrinoCounter
                               << std::endl;

    // Loop over G4 particles
    int particleCounter(0);

    // Find Primary Generator Particles
    std::map<const simb::MCParticle, bool> primaryGeneratorMCParticleMap;
    LArPandoraInput::FindPrimaryParticles(generatorMCParticleVector, primaryGeneratorMCParticleMap);

    for (MCParticleMap::const_iterator iterI = particleMap.begin(), iterEndI = particleMap.end();
         iterI != iterEndI;
         ++iterI) {
      const art::Ptr<simb::MCParticle> particle = iterI->second;

      if (particle->TrackId() != iterI->first)
        throw cet::exception("LArPandora") << "CreatePandoraMCParticles - mc truth information "
                                              "appears to be scrambled in this event";

      if (particle->TrackId() >= settings.m_uidOffset)
        throw cet::exception("LArPandora")
          << "CreatePandoraMCParticles - detected an excessive number of MC particles ("
          << particle->TrackId() << ")";

      ++particleCounter;

      // Find start and end trajectory points
      int firstT(-1), lastT(-1);
      LArPandoraInput::GetTrueStartAndEndPoints(settings, particle, firstT, lastT);

      if (firstT < 0 && lastT < 0) {
        firstT = 0;
        lastT = 0;
      }

      // Lookup position and kinematics at start and end points
      const float vtxX(particle->Vx(firstT));
      const float vtxY(particle->Vy(firstT));
      const float vtxZ(particle->Vz(firstT));

      const float endX(particle->Vx(lastT));
      const float endY(particle->Vy(lastT));
      const float endZ(particle->Vz(lastT));

      const float pX(particle->Px(firstT));
      const float pY(particle->Py(firstT));
      const float pZ(particle->Pz(firstT));
      const float E(particle->E(firstT));

      // Find the source of the mc particle
      int nuanceCode(0);
      const int trackID(particle->TrackId());
      const simb::Origin_t origin(particleInventoryService->TrackIdToMCTruth(trackID).Origin());

      if (LArPandoraInput::IsPrimaryMCParticle(particle, primaryGeneratorMCParticleMap)) {
        nuanceCode = 2001;
      }
      else if (simb::kCosmicRay == origin) {
        nuanceCode = 3000;
      }
      else if (simb::kSingleParticle == origin) {
        nuanceCode = 2000;
      }
      else if (simb::kSuperNovaNeutrino == origin) {
        nuanceCode = 4000;
      }

      // Create 3D Pandora MC Particle
      lar_content::LArMCParticleParameters mcParticleParameters;

      try {
        MCProcessMap processMap;
        FillMCProcessMap(processMap);
        mcParticleParameters.m_nuanceCode = nuanceCode;
        if (processMap.find(particle->Process()) != processMap.end()) {
          mcParticleParameters.m_process = processMap[particle->Process()];
        }
        else {
          mcParticleParameters.m_process = lar_content::MC_PROC_UNKNOWN;
          mf::LogWarning("LArPandora")
            << "CreatePandoraMCParticles - found an unknown process" << std::endl;
        }
        mcParticleParameters.m_energy = E;
        mcParticleParameters.m_particleId = particle->PdgCode();
        mcParticleParameters.m_momentum = pandora::CartesianVector(pX, pY, pZ);
        mcParticleParameters.m_vertex = pandora::CartesianVector(vtxX, vtxY, vtxZ);
        mcParticleParameters.m_endpoint = pandora::CartesianVector(endX, endY, endZ);
        mcParticleParameters.m_mcParticleType = pandora::MC_3D;
        mcParticleParameters.m_pParentAddress = (void*)((intptr_t)particle->TrackId());
      }
      catch (const pandora::StatusCodeException&) {
        mf::LogWarning("LArPandora")
          << "CreatePandoraMCParticles - invalid mc particle parameter provided, all assigned "
             "values must be finite, mc particle omitted "
          << std::endl;
        continue;
      }

      try {
        PANDORA_THROW_RESULT_IF(
          pandora::STATUS_CODE_SUCCESS,
          !=,
          PandoraApi::MCParticle::Create(*pPandora, mcParticleParameters, mcParticleFactory));
      }
      catch (const pandora::StatusCodeException&) {
        mf::LogWarning("LArPandora") << "CreatePandoraMCParticles - unable to create mc particle, "
                                        "insufficient or invalid information supplied "
                                     << std::endl;
        continue;
      }

      // Create Mother/Daughter Links between 3D MC Particles
      const int id_mother(particle->Mother());
      MCParticleMap::const_iterator iterJ = particleMap.find(id_mother);

      if (iterJ != particleMap.end()) {
        try {
          PANDORA_THROW_RESULT_IF(
            pandora::STATUS_CODE_SUCCESS,
            !=,
            PandoraApi::SetMCParentDaughterRelationship(
              *pPandora, (void*)((intptr_t)id_mother), (void*)((intptr_t)particle->TrackId())));
        }
        catch (const pandora::StatusCodeException&) {
          mf::LogWarning("LArPandora") << "CreatePandoraMCParticles - Unable to create mc particle "
                                          "relationship, invalid information supplied "
                                       << std::endl;
          continue;
        }
      }
    }

    mf::LogDebug("LArPandora") << "Number of mc particles: " << particleCounter << std::endl;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraInput::FindPrimaryParticles(
    const RawMCParticleVector& mcParticleVector,
    std::map<const simb::MCParticle, bool>& primaryMCParticleMap)
  {
    for (const simb::MCParticle& mcParticle : mcParticleVector) {
      if ("primary" == mcParticle.Process()) {
        primaryMCParticleMap.emplace(std::make_pair(mcParticle, false));
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  bool LArPandoraInput::IsPrimaryMCParticle(
    const art::Ptr<simb::MCParticle>& mcParticle,
    std::map<const simb::MCParticle, bool>& primaryMCParticleMap)
  {
    for (auto& mcParticleIter : primaryMCParticleMap) {
      if (!mcParticleIter.second) {
        const simb::MCParticle primaryMCParticle(mcParticleIter.first);

        if (std::fabs(primaryMCParticle.Px() - mcParticle->Px()) <
              std::numeric_limits<double>::epsilon() &&
            std::fabs(primaryMCParticle.Py() - mcParticle->Py()) <
              std::numeric_limits<double>::epsilon() &&
            std::fabs(primaryMCParticle.Pz() - mcParticle->Pz()) <
              std::numeric_limits<double>::epsilon()) {
          mcParticleIter.second = true;
          return true;
        }
      }
    }
    return false;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraInput::CreatePandoraMCLinks2D(const Settings& settings,
                                               const IdToHitMap& idToHitMap,
                                               const HitsToTrackIDEs& hitToParticleMap)
  {
    mf::LogDebug("LArPandora") << " *** LArPandoraInput::CreatePandoraMCLinks(...) *** "
                               << std::endl;

    if (!settings.m_pPrimaryPandora)
      throw cet::exception("LArPandora")
        << "CreatePandoraMCLinks2D - primary Pandora instance does not exist ";

    const pandora::Pandora* pPandora(settings.m_pPrimaryPandora);

    for (IdToHitMap::const_iterator iterI = idToHitMap.begin(), iterEndI = idToHitMap.end();
         iterI != iterEndI;
         ++iterI) {
      const int hitID(iterI->first);
      const art::Ptr<recob::Hit> hit(iterI->second);
      //  const geo::WireID hit_WireID(hit->WireID());

      // Get list of associated MC particles
      HitsToTrackIDEs::const_iterator iterJ = hitToParticleMap.find(hit);

      if (hitToParticleMap.end() == iterJ) continue;

      const TrackIDEVector& trackCollection = iterJ->second;

      if (trackCollection.size() == 0)
        throw cet::exception("LArPandora")
          << "CreatePandoraMCLinks2D - found a hit without any associated MC truth information ";

      // Create links between hits and MC particles
      for (unsigned int k = 0; k < trackCollection.size(); ++k) {
        const sim::TrackIDE trackIDE(trackCollection.at(k));
        const int trackID(std::abs(trackIDE.trackID)); // TODO: Find out why std::abs is needed
        const float energyFrac(trackIDE.energyFrac);

        try {
          PANDORA_THROW_RESULT_IF(
            pandora::STATUS_CODE_SUCCESS,
            !=,
            PandoraApi::SetCaloHitToMCParticleRelationship(
              *pPandora, (void*)((intptr_t)hitID), (void*)((intptr_t)trackID), energyFrac));
        }
        catch (const pandora::StatusCodeException&) {
          mf::LogWarning("LArPandora") << "CreatePandoraMCLinks2D - unable to create calo hit to "
                                          "mc particle relationship, invalid information supplied "
                                       << std::endl;
          continue;
        }
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraInput::GetTrueStartAndEndPoints(const Settings& settings,
                                                 const art::Ptr<simb::MCParticle>& particle,
                                                 int& firstT,
                                                 int& lastT)
  {
    art::ServiceHandle<geo::Geometry const> theGeometry;
    firstT = -1;
    lastT = -1;

    for (auto const& tpcid : theGeometry->Iterate<geo::TPCID>()) {
      int thisfirstT(-1), thislastT(-1);
      LArPandoraInput::GetTrueStartAndEndPoints(tpcid, particle, thisfirstT, thislastT);

      if (thisfirstT < 0) continue;

      if (firstT < 0 || thisfirstT < firstT) firstT = thisfirstT;

      if (lastT < 0 || thislastT > lastT) lastT = thislastT;
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraInput::GetTrueStartAndEndPoints(const geo::TPCID& ref_tpcid,
                                                 const art::Ptr<simb::MCParticle>& particle,
                                                 int& startT,
                                                 int& endT)
  {
    art::ServiceHandle<geo::Geometry const> theGeometry;

    bool foundStartPosition(false);
    const int numTrajectoryPoints(static_cast<int>(particle->NumberTrajectoryPoints()));

    for (int nt = 0; nt < numTrajectoryPoints; ++nt) {
      const geo::Point_t pos{particle->Vx(nt), particle->Vy(nt), particle->Vz(nt)};
      geo::TPCID tpcID = theGeometry->FindTPCAtPosition(pos);

      if (!tpcID.isValid) continue;

      if (tpcID != ref_tpcid) continue;

      endT = nt;

      if (!foundStartPosition) {
        startT = endT;
        foundStartPosition = true;
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  float LArPandoraInput::GetTrueX0(const art::Event& e,
                                   const art::Ptr<simb::MCParticle>& particle,
                                   const int nt)
  {
    art::ServiceHandle<geo::Geometry const> theGeometry;
    auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
    auto const det_prop =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clock_data);

    geo::Point_t const pos{particle->Vx(nt), particle->Vy(nt), particle->Vz(nt)};
    auto const tpcID = theGeometry->PositionToTPCID(pos);

    const float vtxT(particle->T(nt));
    const float vtxTDC(clock_data.TPCG4Time2Tick(vtxT));
    const float vtxTDC0(trigger_offset(clock_data));

    auto const [_, sign] = theGeometry->TPC(tpcID).DriftAxisWithSign();
    const float driftDir(sign == geo::DriftSign::Negative ? +1.0 : -1.0);
    return (driftDir * (vtxTDC - vtxTDC0) * det_prop.GetXTicksCoefficient());
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  double LArPandoraInput::GetMips(detinfo::DetectorPropertiesData const& detProp,
                                  const Settings& settings,
                                  const double hit_Charge,
                                  const geo::View_t hit_View)
  {
    auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout const>()->Get();

    // TODO: Unite this procedure with other calorimetry procedures under development
    const double dQdX(hit_Charge / (wireReadoutGeom.Plane({0, 0}, hit_View).WirePitch())); // ADC/cm
    const double dQdX_e(dQdX /
                        (detProp.ElectronsToADC() * settings.m_recombination_factor)); // e/cm
    const double dEdX(settings.m_useBirksCorrection ?
                        detProp.BirksCorrection(dQdX_e) :
                        dQdX_e * 1000. / util::kGeVToElectrons); // MeV/cm
    double mips(dEdX / settings.m_dEdX_mip);

    if (mips < 0.) mips = settings.m_mips_if_negative;

    if (mips > settings.m_mips_max) mips = settings.m_mips_max;

    return mips;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraInput::FillMCProcessMap(MCProcessMap& processMap)
  {
    // QGSP_BERT and EM standard physics list mappings
    processMap["unknown"] = lar_content::MC_PROC_UNKNOWN;
    processMap["primary"] = lar_content::MC_PROC_PRIMARY;
    processMap["compt"] = lar_content::MC_PROC_COMPT;
    processMap["phot"] = lar_content::MC_PROC_PHOT;
    processMap["annihil"] = lar_content::MC_PROC_ANNIHIL;
    processMap["eIoni"] = lar_content::MC_PROC_E_IONI;
    processMap["eBrem"] = lar_content::MC_PROC_E_BREM;
    processMap["conv"] = lar_content::MC_PROC_CONV;
    processMap["muIoni"] = lar_content::MC_PROC_MU_IONI;
    processMap["muMinusCaptureAtRest"] = lar_content::MC_PROC_MU_MINUS_CAPTURE_AT_REST;
    processMap["neutronInelastic"] = lar_content::MC_PROC_NEUTRON_INELASTIC;
    processMap["nCapture"] = lar_content::MC_PROC_N_CAPTURE;
    processMap["hadElastic"] = lar_content::MC_PROC_HAD_ELASTIC;
    processMap["Decay"] = lar_content::MC_PROC_DECAY;
    processMap["CoulombScat"] = lar_content::MC_PROC_COULOMB_SCAT;
    processMap["muBrems"] = lar_content::MC_PROC_MU_BREM;
    processMap["muPairProd"] = lar_content::MC_PROC_MU_PAIR_PROD;
    processMap["PhotonInelastic"] = lar_content::MC_PROC_PHOTON_INELASTIC;
    processMap["hIoni"] = lar_content::MC_PROC_HAD_IONI;
    processMap["protonInelastic"] = lar_content::MC_PROC_PROTON_INELASTIC;
    processMap["pi+Inelastic"] = lar_content::MC_PROC_PI_PLUS_INELASTIC;
    processMap["CHIPSNuclearCaptureAtRest"] = lar_content::MC_PROC_CHIPS_NUCLEAR_CAPTURE_AT_REST;
    processMap["pi-Inelastic"] = lar_content::MC_PROC_PI_MINUS_INELASTIC;
    processMap["Transportation"] = lar_content::MC_PROC_TRANSPORTATION;
    processMap["Rayl"] = lar_content::MC_PROC_RAYLEIGH;
    processMap["hBrems"] = lar_content::MC_PROC_HAD_BREM;
    processMap["hPairProd"] = lar_content::MC_PROC_HAD_PAIR_PROD;
    processMap["ionIoni"] = lar_content::MC_PROC_ION_IONI;
    processMap["nKiller"] = lar_content::MC_PROC_NEUTRON_KILLER;
    processMap["ionInelastic"] = lar_content::MC_PROC_ION_INELASTIC;
    processMap["He3Inelastic"] = lar_content::MC_PROC_HE3_INELASTIC;
    processMap["alphaInelastic"] = lar_content::MC_PROC_ALPHA_INELASTIC;
    processMap["anti_He3Inelastic"] = lar_content::MC_PROC_ANTI_HE3_INELASTIC;
    processMap["anti_alphaInelastic"] = lar_content::MC_PROC_ANTI_ALPHA_INELASTIC;
    processMap["hFritiofCaptureAtRest"] = lar_content::MC_PROC_HAD_FRITIOF_CAPTURE_AT_REST;
    processMap["anti_deuteronInelastic"] = lar_content::MC_PROC_ANTI_DEUTERON_INELASTIC;
    processMap["anti_neutronInelastic"] = lar_content::MC_PROC_ANTI_NEUTRON_INELASTIC;
    processMap["anti_protonInelastic"] = lar_content::MC_PROC_ANTI_PROTON_INELASTIC;
    processMap["anti_tritonInelastic"] = lar_content::MC_PROC_ANTI_TRITON_INELASTIC;
    processMap["dInelastic"] = lar_content::MC_PROC_DEUTERON_INELASTIC;
    processMap["electronNuclear"] = lar_content::MC_PROC_ELECTRON_NUCLEAR;
    processMap["photonNuclear"] = lar_content::MC_PROC_PHOTON_NUCLEAR;
    processMap["kaon+Inelastic"] = lar_content::MC_PROC_KAON_PLUS_INELASTIC;
    processMap["kaon-Inelastic"] = lar_content::MC_PROC_KAON_MINUS_INELASTIC;
    processMap["hBertiniCaptureAtRest"] = lar_content::MC_PROC_HAD_BERTINI_CAPTURE_AT_REST;
    processMap["lambdaInelastic"] = lar_content::MC_PROC_LAMBDA_INELASTIC;
    processMap["muonNuclear"] = lar_content::MC_PROC_MU_NUCLEAR;
    processMap["tInelastic"] = lar_content::MC_PROC_TRITON_INELASTIC;
    processMap["primaryBackground"] = lar_content::MC_PROC_PRIMARY_BACKGROUND;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------------------------------------------------

  LArPandoraInput::Settings::Settings()
    : m_pPrimaryPandora(nullptr)
    , m_useHitWidths(true)
    , m_useBirksCorrection(false)
    , m_useActiveBoundingBox(false)
    , m_uidOffset(100000000)
    , m_hitCounterOffset(0)
    , m_dx_cm(0.5)
    , m_int_cm(84.)
    , m_rad_cm(14.)
    , m_dEdX_mip(2.)
    , m_mips_max(50.)
    , m_mips_if_negative(0.)
    , m_mips_to_gev(3.5e-4)
    , m_recombination_factor(0.63)
  {}

} // namespace lar_pandora
