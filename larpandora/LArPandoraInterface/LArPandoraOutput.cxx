/**
 *  @file   larpandora/LArPandoraInterface/LArPandoraOutput.cxx
 *
 *  @brief  Helper functions for processing outputs from pandora
 *
 */
#include "larpandora/LArPandoraInterface/LArPandoraOutput.h"

#include "art/Framework/Principal/Event.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "larreco/ClusterFinder/ClusterCreator.h"
#include "larreco/RecoAlg/ClusterParamsImportWrapper.h"
#include "larreco/RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"

#include "Objects/CaloHit.h"
#include "Objects/Cluster.h"
#include "Objects/ParticleFlowObject.h"
#include "Objects/Vertex.h"

#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric> // std::iota()

namespace lar_pandora {

  void LArPandoraOutput::ProduceArtOutput(const Settings& settings,
                                          const IdToHitMap& idToHitMap,
                                          art::Event& evt)
  {
    settings.Validate();
    const std::string instanceLabel(
      settings.m_shouldProduceAllOutcomes ? settings.m_allOutcomesInstanceLabel : "");
    const std::string testBeamInteractionVertexInstanceLabel(
      instanceLabel + settings.m_testBeamInteractionVerticesInstanceLabel);

    // Set up mandatory output collections
    PFParticleCollection outputParticles(new std::vector<recob::PFParticle>);
    VertexCollection outputVertices(new std::vector<recob::Vertex>);
    ClusterCollection outputClusters(new std::vector<recob::Cluster>);
    SpacePointCollection outputSpacePoints(new std::vector<recob::SpacePoint>);
    PFParticleMetadataCollection outputParticleMetadata(
      new std::vector<larpandoraobj::PFParticleMetadata>);

    // Set up optional output collections
    VertexCollection outputTestBeamInteractionVertices(
      settings.m_shouldProduceTestBeamInteractionVertices ? new std::vector<recob::Vertex> :
                                                            nullptr);
    T0Collection outputT0s(settings.m_shouldRunStitching ? new std::vector<anab::T0> : nullptr);
    SliceCollection outputSlices(settings.m_shouldProduceSlices ? new std::vector<recob::Slice> :
                                                                  nullptr);

    // Set up mandatory output associations
    PFParticleToMetadataCollection outputParticlesToMetadata(
      new art::Assns<recob::PFParticle, larpandoraobj::PFParticleMetadata>);
    PFParticleToSpacePointCollection outputParticlesToSpacePoints(
      new art::Assns<recob::PFParticle, recob::SpacePoint>);
    PFParticleToClusterCollection outputParticlesToClusters(
      new art::Assns<recob::PFParticle, recob::Cluster>);
    PFParticleToVertexCollection outputParticlesToVertices(
      new art::Assns<recob::PFParticle, recob::Vertex>);
    ClusterToHitCollection outputClustersToHits(new art::Assns<recob::Cluster, recob::Hit>);
    SpacePointToHitCollection outputSpacePointsToHits(
      new art::Assns<recob::SpacePoint, recob::Hit>);
    SliceToHitCollection outputSlicesToHits(new art::Assns<recob::Slice, recob::Hit>);

    // Set up optional output associations
    PFParticleToVertexCollection outputParticlesToTestBeamInteractionVertices(
      settings.m_shouldProduceTestBeamInteractionVertices ?
        new art::Assns<recob::PFParticle, recob::Vertex> :
        nullptr);
    PFParticleToT0Collection outputParticlesToT0s(
      settings.m_shouldRunStitching ? new art::Assns<recob::PFParticle, anab::T0> : nullptr);
    PFParticleToSliceCollection outputParticlesToSlices(
      settings.m_shouldProduceSlices ? new art::Assns<recob::PFParticle, recob::Slice> : nullptr);

    // Collect immutable lists of pandora collections that we should convert to ART format
    const pandora::PfoVector pfoVector(
      settings.m_shouldProduceAllOutcomes ?
        LArPandoraOutput::CollectAllPfoOutcomes(settings.m_pPrimaryPandora) :
        LArPandoraOutput::CollectPfos(settings.m_pPrimaryPandora));

    IdToIdVectorMap pfoToVerticesMap, pfoToTestBeamInteractionVerticesMap;
    const pandora::VertexVector vertexVector(LArPandoraOutput::CollectVertices(
      pfoVector, pfoToVerticesMap, lar_content::LArPfoHelper::GetVertex));
    const pandora::VertexVector testBeamInteractionVertexVector(
      settings.m_shouldProduceTestBeamInteractionVertices ?
        LArPandoraOutput::CollectVertices(pfoVector,
                                          pfoToTestBeamInteractionVerticesMap,
                                          lar_content::LArPfoHelper::GetTestBeamInteractionVertex) :
        pandora::VertexVector());

    IdToIdVectorMap pfoToClustersMap;
    const pandora::ClusterList clusterList(
      LArPandoraOutput::CollectClusters(pfoVector, pfoToClustersMap));

    IdToIdVectorMap pfoToThreeDHitsMap;
    const pandora::CaloHitList threeDHitList(
      LArPandoraOutput::Collect3DHits(pfoVector, pfoToThreeDHitsMap));

    // Get mapping from pandora hits to art hits
    CaloHitToArtHitMap pandoraHitToArtHitMap;
    LArPandoraOutput::GetPandoraToArtHitMap(
      clusterList, threeDHitList, idToHitMap, pandoraHitToArtHitMap);

    // Build the ART outputs from the pandora objects
    LArPandoraOutput::BuildVertices(vertexVector, outputVertices);

    if (settings.m_shouldProduceTestBeamInteractionVertices)
      LArPandoraOutput::BuildVertices(testBeamInteractionVertexVector,
                                      outputTestBeamInteractionVertices);

    LArPandoraOutput::BuildSpacePoints(evt,
                                       instanceLabel,
                                       threeDHitList,
                                       pandoraHitToArtHitMap,
                                       outputSpacePoints,
                                       outputSpacePointsToHits);

    IdToIdVectorMap pfoToArtClustersMap;
    LArPandoraOutput::BuildClusters(evt,
                                    instanceLabel,
                                    clusterList,
                                    pandoraHitToArtHitMap,
                                    pfoToClustersMap,
                                    outputClusters,
                                    outputClustersToHits,
                                    pfoToArtClustersMap);

    LArPandoraOutput::BuildPFParticles(evt,
                                       instanceLabel,
                                       pfoVector,
                                       pfoToVerticesMap,
                                       pfoToThreeDHitsMap,
                                       pfoToArtClustersMap,
                                       outputParticles,
                                       outputParticlesToVertices,
                                       outputParticlesToSpacePoints,
                                       outputParticlesToClusters);

    LArPandoraOutput::BuildParticleMetadata(
      evt, instanceLabel, pfoVector, outputParticleMetadata, outputParticlesToMetadata);

    if (settings.m_shouldProduceSlices)
      LArPandoraOutput::BuildSlices(settings,
                                    settings.m_pPrimaryPandora,
                                    evt,
                                    instanceLabel,
                                    pfoVector,
                                    idToHitMap,
                                    outputSlices,
                                    outputParticlesToSlices,
                                    outputSlicesToHits);

    if (settings.m_shouldRunStitching)
      LArPandoraOutput::BuildT0s(evt, instanceLabel, pfoVector, outputT0s, outputParticlesToT0s);

    if (settings.m_shouldProduceTestBeamInteractionVertices)
      LArPandoraOutput::AssociateAdditionalVertices(evt,
                                                    instanceLabel,
                                                    pfoVector,
                                                    pfoToTestBeamInteractionVerticesMap,
                                                    outputParticlesToTestBeamInteractionVertices);

    // Add the outputs to the event
    evt.put(std::move(outputParticles), instanceLabel);
    evt.put(std::move(outputSpacePoints), instanceLabel);
    evt.put(std::move(outputClusters), instanceLabel);
    evt.put(std::move(outputVertices), instanceLabel);
    evt.put(std::move(outputParticleMetadata), instanceLabel);

    evt.put(std::move(outputParticlesToMetadata), instanceLabel);
    evt.put(std::move(outputParticlesToSpacePoints), instanceLabel);
    evt.put(std::move(outputParticlesToClusters), instanceLabel);
    evt.put(std::move(outputParticlesToVertices), instanceLabel);
    evt.put(std::move(outputParticlesToSlices), instanceLabel);
    evt.put(std::move(outputSpacePointsToHits), instanceLabel);
    evt.put(std::move(outputClustersToHits), instanceLabel);

    if (settings.m_shouldProduceTestBeamInteractionVertices) {
      evt.put(std::move(outputTestBeamInteractionVertices), testBeamInteractionVertexInstanceLabel);
      evt.put(std::move(outputParticlesToTestBeamInteractionVertices),
              testBeamInteractionVertexInstanceLabel);
    }

    if (settings.m_shouldRunStitching) {
      evt.put(std::move(outputT0s), instanceLabel);
      evt.put(std::move(outputParticlesToT0s), instanceLabel);
    }

    if (settings.m_shouldProduceSlices) {
      evt.put(std::move(outputSlices), instanceLabel);
      evt.put(std::move(outputSlicesToHits), instanceLabel);
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  bool LArPandoraOutput::GetPandoraInstance(const pandora::Pandora* const pPrimaryPandora,
                                            const std::string& name,
                                            const pandora::Pandora*& pPandoraInstance)
  {
    if (!pPrimaryPandora)
      throw cet::exception("LArPandora") << " LArPandoraOutput::GetPandoraInstance--- input "
                                            "primary pandora instance address is invalid ";

    if (pPandoraInstance)
      throw cet::exception("LArPandora") << " LArPandoraOutput::GetPandoraInstance--- the input "
                                            "pandora instance address is non-null ";

    for (const pandora::Pandora* const pPandora :
         MultiPandoraApi::GetDaughterPandoraInstanceList(pPrimaryPandora)) {
      if (pPandora->GetName() != name) continue;

      if (pPandoraInstance)
        throw cet::exception("LArPandora")
          << " LArPandoraOutput::GetPandoraInstance--- found multiple pandora instances with name: "
          << name;

      pPandoraInstance = pPandora;
    }

    return static_cast<bool>(pPandoraInstance);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraOutput::GetPandoraSlices(const pandora::Pandora* const pPrimaryPandora,
                                          pandora::PfoVector& slicePfos)
  {
    if (!slicePfos.empty())
      throw cet::exception("LArPandora")
        << " LArPandoraOutput::GetPandoraSlices--- Input slicePfo vector is not empty ";

    // Get the pandora slicing worker - if it exists
    const pandora::Pandora* pSlicingWorker(nullptr);
    if (!LArPandoraOutput::GetPandoraInstance(pPrimaryPandora, "SlicingWorker", pSlicingWorker))
      return;

    // Get the slice PFOs - one PFO per slice
    const pandora::PfoList* pSlicePfoList(nullptr);
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS,
                            !=,
                            PandoraApi::GetCurrentPfoList(*pSlicingWorker, pSlicePfoList));

    slicePfos.insert(slicePfos.end(), pSlicePfoList->begin(), pSlicePfoList->end());
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  pandora::PfoVector LArPandoraOutput::CollectAllPfoOutcomes(
    const pandora::Pandora* const pPrimaryPandora)
  {
    pandora::PfoList collectedPfos;

    // Get the list of slice pfos - one per slice
    pandora::PfoVector slicePfos;
    LArPandoraOutput::GetPandoraSlices(pPrimaryPandora, slicePfos);

    // Identify the pandora worker instances by their name
    const pandora::Pandora* pSliceNuWorker(nullptr);
    if (!LArPandoraOutput::GetPandoraInstance(pPrimaryPandora, "SliceNuWorker", pSliceNuWorker))
      throw cet::exception("LArPandora")
        << " LArPandoraOutput::CollectAllPfoOutcomes--- Can't find slice nu worker instance. ";

    const pandora::Pandora* pSliceCRWorker(nullptr);
    if (!LArPandoraOutput::GetPandoraInstance(pPrimaryPandora, "SliceCRWorker", pSliceCRWorker))
      throw cet::exception("LArPandora")
        << " LArPandoraOutput::CollectAllPfoOutcomes--- Can't find slice CR worker instance. ";

    // Collect slices under both reconstruction hypotheses
    for (unsigned int sliceIndex = 0; sliceIndex < slicePfos.size(); ++sliceIndex) {
      const pandora::PfoList* pNuPfoList(nullptr);
      if (pandora::STATUS_CODE_SUCCESS ==
          PandoraApi::GetPfoList(
            *pSliceNuWorker, "NeutrinoParticles3D" + std::to_string(sliceIndex), pNuPfoList))
        collectedPfos.insert(collectedPfos.end(), pNuPfoList->begin(), pNuPfoList->end());

      const pandora::PfoList* pCRPfoList(nullptr);
      if (pandora::STATUS_CODE_SUCCESS ==
          PandoraApi::GetPfoList(
            *pSliceCRWorker, "MuonParticles3D" + std::to_string(sliceIndex), pCRPfoList))
        collectedPfos.insert(collectedPfos.end(), pCRPfoList->begin(), pCRPfoList->end());
    }

    // Get the list of the parent pfos from the primary pandora instance
    const pandora::PfoList* pParentPfoList(nullptr);
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS,
                            !=,
                            PandoraApi::GetCurrentPfoList(*pPrimaryPandora, pParentPfoList));

    // Collect clear cosmic-rays
    for (const pandora::ParticleFlowObject* const pPfo : *pParentPfoList) {
      if (LArPandoraOutput::IsClearCosmic(pPfo)) collectedPfos.push_back(pPfo);
    }

    // Collect all pfos that are downstream of the parents we have collected
    pandora::PfoVector pfoVector;
    LArPandoraOutput::CollectPfos(collectedPfos, pfoVector);

    return pfoVector;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  bool LArPandoraOutput::IsClearCosmic(const pandora::ParticleFlowObject* const pPfo)
  {
    const pandora::ParticleFlowObject* const pParent(lar_content::LArPfoHelper::GetParentPfo(pPfo));

    const auto& properties(pParent->GetPropertiesMap());
    const auto it(properties.find("IsClearCosmic"));

    if (it == properties.end()) return false;

    return static_cast<bool>(std::round(it->second));
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  bool LArPandoraOutput::IsFromSlice(const pandora::ParticleFlowObject* const pPfo)
  {
    const pandora::ParticleFlowObject* const pParent(lar_content::LArPfoHelper::GetParentPfo(pPfo));

    const auto& properties(pParent->GetPropertiesMap());
    return (properties.find("SliceIndex") != properties.end());
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  unsigned int LArPandoraOutput::GetSliceIndex(const pandora::ParticleFlowObject* const pPfo)
  {
    const pandora::ParticleFlowObject* const pParent(lar_content::LArPfoHelper::GetParentPfo(pPfo));

    const auto& properties(pParent->GetPropertiesMap());
    const auto it(properties.find("SliceIndex"));

    if (it == properties.end())
      throw cet::exception("LArPandora")
        << " LArPandoraOutput::GetSliceIndex--- Input PFO was not from a slice ";

    return static_cast<unsigned int>(std::round(it->second));
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  pandora::PfoVector LArPandoraOutput::CollectPfos(const pandora::Pandora* const pPrimaryPandora)
  {
    const pandora::PfoList* pParentPfoList(nullptr);
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS,
                            !=,
                            PandoraApi::GetCurrentPfoList(*pPrimaryPandora, pParentPfoList));

    pandora::PfoVector pfoVector;
    LArPandoraOutput::CollectPfos(*pParentPfoList, pfoVector);

    return pfoVector;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraOutput::CollectPfos(const pandora::PfoList& parentPfoList,
                                     pandora::PfoVector& pfoVector)
  {
    if (!pfoVector.empty())
      throw cet::exception("LArPandora")
        << " LArPandoraOutput::CollectPfos--- trying to collect pfos into a non-empty list ";

    pandora::PfoList pfoList;
    lar_content::LArPfoHelper::GetAllConnectedPfos(parentPfoList, pfoList);

    pfoVector.insert(pfoVector.end(), pfoList.begin(), pfoList.end());
    std::sort(pfoVector.begin(), pfoVector.end(), lar_content::LArPfoHelper::SortByNHits);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  pandora::VertexVector LArPandoraOutput::CollectVertices(
    const pandora::PfoVector& pfoVector,
    IdToIdVectorMap& pfoToVerticesMap,
    std::function<const pandora::Vertex* const(const pandora::ParticleFlowObject* const)> fCriteria)
  {
    pandora::VertexVector vertexVector;

    for (unsigned int pfoId = 0; pfoId < pfoVector.size(); ++pfoId) {
      const pandora::ParticleFlowObject* const pPfo(pfoVector.at(pfoId));

      if (pPfo->GetVertexList().empty()) continue;

      try {
        const pandora::Vertex* const pVertex(fCriteria(pPfo));

        // Get the vertex ID and add it to the vertex list if required
        const auto it(std::find(vertexVector.begin(), vertexVector.end(), pVertex));
        const bool isInList(it != vertexVector.end());
        const size_t vertexId(isInList ? std::distance(vertexVector.begin(), it) :
                                         vertexVector.size());

        if (!isInList) vertexVector.push_back(pVertex);

        if (!pfoToVerticesMap.insert(IdToIdVectorMap::value_type(pfoId, {vertexId})).second)
          throw cet::exception("LArPandora")
            << " LArPandoraOutput::CollectVertices --- repeated pfos in input list ";
      }
      catch (const pandora::StatusCodeException&) {
        continue;
      }
    }

    return vertexVector;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  pandora::ClusterList LArPandoraOutput::CollectClusters(const pandora::PfoVector& pfoVector,
                                                         IdToIdVectorMap& pfoToClustersMap)
  {
    pandora::ClusterList clusterList;

    for (unsigned int pfoId = 0; pfoId < pfoVector.size(); ++pfoId) {
      const pandora::ParticleFlowObject* const pPfo(pfoVector.at(pfoId));

      // Get the sorted list of clusters from the pfo
      pandora::ClusterList clusters;
      lar_content::LArPfoHelper::GetTwoDClusterList(pPfo, clusters);
      clusters.sort(lar_content::LArClusterHelper::SortByNHits);

      // Get incrementing id's for each cluster
      IdVector clusterIds(clusters.size());
      std::iota(clusterIds.begin(), clusterIds.end(), clusterList.size());

      clusterList.insert(clusterList.end(), clusters.begin(), clusters.end());

      if (!pfoToClustersMap.insert(IdToIdVectorMap::value_type(pfoId, clusterIds)).second)
        throw cet::exception("LArPandora")
          << " LArPandoraOutput::CollectClusters --- repeated pfos in input list ";
    }

    return clusterList;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  pandora::CaloHitList LArPandoraOutput::Collect3DHits(const pandora::PfoVector& pfoVector,
                                                       IdToIdVectorMap& pfoToThreeDHitsMap)
  {
    pandora::CaloHitList caloHitList;

    for (unsigned int pfoId = 0; pfoId < pfoVector.size(); ++pfoId) {
      const pandora::ParticleFlowObject* const pPfo(pfoVector.at(pfoId));

      if (!pfoToThreeDHitsMap.insert(IdToIdVectorMap::value_type(pfoId, {})).second)
        throw cet::exception("LArPandora")
          << " LArPandoraOutput::Collect3DHits --- repeated pfos in input list ";

      pandora::CaloHitVector sorted3DHits;
      LArPandoraOutput::Collect3DHits(pPfo, sorted3DHits);

      for (const pandora::CaloHit* const pCaloHit3D : sorted3DHits) {
        if (pandora::TPC_3D !=
            pCaloHit3D
              ->GetHitType()) // TODO decide if this is required, or should I just insert them?
          throw cet::exception("LArPandora")
            << " LArPandoraOutput::Collect3DHits --- found a 2D hit in a 3D cluster";

        pfoToThreeDHitsMap.at(pfoId).push_back(caloHitList.size());
        caloHitList.push_back(pCaloHit3D);
      }
    }

    return caloHitList;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraOutput::Collect3DHits(const pandora::ParticleFlowObject* const pPfo,
                                       pandora::CaloHitVector& caloHits)
  {
    // Get the sorted list of 3D hits associated with the pfo
    pandora::CaloHitList threeDHits;
    lar_content::LArPfoHelper::GetCaloHits(pPfo, pandora::TPC_3D, threeDHits);

    caloHits.insert(caloHits.end(), threeDHits.begin(), threeDHits.end());
    std::sort(caloHits.begin(), caloHits.end(), lar_content::LArClusterHelper::SortHitsByPosition);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraOutput::GetPandoraToArtHitMap(const pandora::ClusterList& clusterList,
                                               const pandora::CaloHitList& threeDHitList,
                                               const IdToHitMap& idToHitMap,
                                               CaloHitToArtHitMap& pandoraHitToArtHitMap)
  {
    // Collect 2D hits from clusters
    for (const pandora::Cluster* const pCluster : clusterList) {
      if (pandora::TPC_3D == lar_content::LArClusterHelper::GetClusterHitType(pCluster))
        throw cet::exception("LArPandora")
          << " LArPandoraOutput::GetPandoraToArtHitMap --- found a 3D input cluster ";

      pandora::CaloHitVector sortedHits;
      LArPandoraOutput::GetHitsInCluster(pCluster, sortedHits);

      for (const pandora::CaloHit* const pCaloHit : sortedHits) {
        if (!pandoraHitToArtHitMap
               .insert(CaloHitToArtHitMap::value_type(
                 pCaloHit, LArPandoraOutput::GetHit(idToHitMap, pCaloHit)))
               .second)
          throw cet::exception("LArPandora")
            << " LArPandoraOutput::GetPandoraToArtHitMap --- found repeated input hits ";
      }
    }

    for (const pandora::CaloHit* const pCaloHit : threeDHitList) {
      if (pCaloHit->GetHitType() != pandora::TPC_3D)
        throw cet::exception("LArPandora")
          << " LArPandoraOutput::GetPandoraToArtHitMap --- found a non-3D hit in the input list ";

      // ATTN get the 2D calo hit from the 3D calo hit then find the art hit!
      if (!pandoraHitToArtHitMap
             .insert(CaloHitToArtHitMap::value_type(
               pCaloHit,
               LArPandoraOutput::GetHit(
                 idToHitMap, static_cast<const pandora::CaloHit*>(pCaloHit->GetParentAddress()))))
             .second)
        throw cet::exception("LArPandora")
          << " LArPandoraOutput::GetPandoraToArtHitMap --- found repeated input hits ";
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  art::Ptr<recob::Hit> LArPandoraOutput::GetHit(const IdToHitMap& idToHitMap,
                                                const pandora::CaloHit* const pCaloHit)
  {
    //  TODO make this less evil

    // ATTN The CaloHit can come from the primary pandora instance (depth = 0) or one of its daughers (depth = 1).
    //      Here we keep trying to access the ART hit increasing the depth step-by-step
    for (unsigned int depth = 0, maxDepth = 2; depth < maxDepth; ++depth) {
      // Navigate to the hit address in the pandora master instance (assuming the depth is correct)
      const pandora::CaloHit* pParentCaloHit = pCaloHit;
      for (unsigned int i = 0; i < depth; ++i)
        pParentCaloHit = static_cast<const pandora::CaloHit*>(pCaloHit->GetParentAddress());

      // Attempt to find the mapping from the "parent" calo hit to the ART hit
      const void* const pHitAddress(pParentCaloHit->GetParentAddress());
      const intptr_t hitID_temp((intptr_t)(pHitAddress));
      const int hitID((int)(hitID_temp));

      IdToHitMap::const_iterator artIter = idToHitMap.find(hitID);

      // If there is no such mapping from "parent" calo hit to the ART hit, then increase the depth and try again!
      if (idToHitMap.end() == artIter) continue;

      return artIter->second;
    }

    throw cet::exception("LArPandora")
      << " LArPandoraOutput::GetHit --- found a Pandora hit without a parent ART hit ";
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraOutput::BuildVertices(const pandora::VertexVector& vertexVector,
                                       VertexCollection& outputVertices)
  {
    for (size_t vertexId = 0; vertexId < vertexVector.size(); ++vertexId)
      outputVertices->push_back(LArPandoraOutput::BuildVertex(vertexVector.at(vertexId), vertexId));
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraOutput::BuildSpacePoints(const art::Event& event,
                                          const std::string& instanceLabel,
                                          const pandora::CaloHitList& threeDHitList,
                                          const CaloHitToArtHitMap& pandoraHitToArtHitMap,
                                          SpacePointCollection& outputSpacePoints,
                                          SpacePointToHitCollection& outputSpacePointsToHits)
  {
    pandora::CaloHitVector threeDHitVector;
    threeDHitVector.insert(threeDHitVector.end(), threeDHitList.begin(), threeDHitList.end());

    for (unsigned int hitId = 0; hitId < threeDHitVector.size(); hitId++) {
      const pandora::CaloHit* const pCaloHit(threeDHitVector.at(hitId));

      CaloHitToArtHitMap::const_iterator it(pandoraHitToArtHitMap.find(pCaloHit));
      if (it == pandoraHitToArtHitMap.end())
        throw cet::exception("LArPandora") << " LArPandoraOutput::BuildSpacePoints --- found a "
                                              "pandora hit without a corresponding art hit ";

      LArPandoraOutput::AddAssociation(
        event, instanceLabel, hitId, {it->second}, outputSpacePointsToHits);
      outputSpacePoints->push_back(LArPandoraOutput::BuildSpacePoint(pCaloHit, hitId));
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraOutput::BuildClusters(const art::Event& event,
                                       const std::string& instanceLabel,
                                       const pandora::ClusterList& clusterList,
                                       const CaloHitToArtHitMap& pandoraHitToArtHitMap,
                                       const IdToIdVectorMap& pfoToClustersMap,
                                       ClusterCollection& outputClusters,
                                       ClusterToHitCollection& outputClustersToHits,
                                       IdToIdVectorMap& pfoToArtClustersMap)
  {
    cluster::StandardClusterParamsAlg clusterParamAlgo;

    art::ServiceHandle<geo::Geometry const> geom{};
    auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout const>()->Get();
    auto const clock_data =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
    auto const det_prop =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clock_data);
    util::GeometryUtilities const gser{*geom, wireReadoutGeom, clock_data, det_prop};

    // Produce the art clusters
    size_t nextClusterId(0);
    IdToIdVectorMap pandoraClusterToArtClustersMap;
    for (const pandora::Cluster* const pCluster : clusterList) {
      std::vector<HitVector> hitVectors;
      const std::vector<recob::Cluster> clusters(
        LArPandoraOutput::BuildClusters(gser,
                                        pCluster,
                                        clusterList,
                                        pandoraHitToArtHitMap,
                                        pandoraClusterToArtClustersMap,
                                        hitVectors,
                                        nextClusterId,
                                        clusterParamAlgo));

      if (hitVectors.size() != clusters.size())
        throw cet::exception("LArPandora")
          << " LArPandoraOutput::BuildClusters --- invalid hit vectors for clusters produced ";

      for (unsigned int i = 0; i < clusters.size(); ++i) {
        LArPandoraOutput::AddAssociation(
          event, instanceLabel, nextClusterId - 1, hitVectors.at(i), outputClustersToHits);
        outputClusters->push_back(clusters.at(i));
      }
    }

    // Get mapping from pfo id to art cluster id
    for (IdToIdVectorMap::const_iterator it = pfoToClustersMap.begin();
         it != pfoToClustersMap.end();
         ++it) {
      if (!pfoToArtClustersMap.insert(IdToIdVectorMap::value_type(it->first, {})).second)
        throw cet::exception("LArPandora")
          << " LArPandoraOutput::BuildClusters --- repeated pfo ids ";

      for (const size_t pandoraClusterId : it->second) {
        IdToIdVectorMap::const_iterator it2(pandoraClusterToArtClustersMap.find(pandoraClusterId));

        if (it2 == pandoraClusterToArtClustersMap.end())
          throw cet::exception("LArPandora") << " LArPandoraOutput::BuildClusters --- found a "
                                                "pandora cluster with no associated recob cluster ";

        for (const size_t recobClusterId : it2->second)
          pfoToArtClustersMap.at(it->first).push_back(recobClusterId);
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraOutput::BuildPFParticles(
    const art::Event& event,
    const std::string& instanceLabel,
    const pandora::PfoVector& pfoVector,
    const IdToIdVectorMap& pfoToVerticesMap,
    const IdToIdVectorMap& pfoToThreeDHitsMap,
    const IdToIdVectorMap& pfoToArtClustersMap,
    PFParticleCollection& outputParticles,
    PFParticleToVertexCollection& outputParticlesToVertices,
    PFParticleToSpacePointCollection& outputParticlesToSpacePoints,
    PFParticleToClusterCollection& outputParticlesToClusters)
  {
    for (unsigned int pfoId = 0; pfoId < pfoVector.size(); ++pfoId) {
      const pandora::ParticleFlowObject* const pPfo(pfoVector.at(pfoId));

      outputParticles->push_back(LArPandoraOutput::BuildPFParticle(pPfo, pfoId, pfoVector));

      // Associations from PFParticle
      if (pfoToVerticesMap.find(pfoId) != pfoToVerticesMap.end())
        LArPandoraOutput::AddAssociation(
          event, instanceLabel, pfoId, pfoToVerticesMap, outputParticlesToVertices);

      if (pfoToThreeDHitsMap.find(pfoId) != pfoToThreeDHitsMap.end())
        LArPandoraOutput::AddAssociation(
          event, instanceLabel, pfoId, pfoToThreeDHitsMap, outputParticlesToSpacePoints);

      if (pfoToArtClustersMap.find(pfoId) != pfoToArtClustersMap.end())
        LArPandoraOutput::AddAssociation(
          event, instanceLabel, pfoId, pfoToArtClustersMap, outputParticlesToClusters);
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraOutput::AssociateAdditionalVertices(
    const art::Event& event,
    const std::string& instanceLabel,
    const pandora::PfoVector& pfoVector,
    const IdToIdVectorMap& pfoToVerticesMap,
    PFParticleToVertexCollection& outputParticlesToVertices)
  {
    for (unsigned int pfoId = 0; pfoId < pfoVector.size(); ++pfoId) {
      if (pfoToVerticesMap.find(pfoId) != pfoToVerticesMap.end())
        LArPandoraOutput::AddAssociation(
          event, instanceLabel, pfoId, pfoToVerticesMap, outputParticlesToVertices);
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraOutput::BuildParticleMetadata(
    const art::Event& event,
    const std::string& instanceLabel,
    const pandora::PfoVector& pfoVector,
    PFParticleMetadataCollection& outputParticleMetadata,
    PFParticleToMetadataCollection& outputParticlesToMetadata)
  {
    for (unsigned int pfoId = 0; pfoId < pfoVector.size(); ++pfoId) {
      const pandora::ParticleFlowObject* const pPfo(pfoVector.at(pfoId));

      LArPandoraOutput::AddAssociation(
        event, instanceLabel, pfoId, outputParticleMetadata->size(), outputParticlesToMetadata);
      larpandoraobj::PFParticleMetadata pPFParticleMetadata(
        LArPandoraHelper::GetPFParticleMetadata(pPfo));
      outputParticleMetadata->push_back(pPFParticleMetadata);
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraOutput::BuildSlices(const Settings& settings,
                                     const pandora::Pandora* const pPrimaryPandora,
                                     const art::Event& event,
                                     const std::string& instanceLabel,
                                     const pandora::PfoVector& pfoVector,
                                     const IdToHitMap& idToHitMap,
                                     SliceCollection& outputSlices,
                                     PFParticleToSliceCollection& outputParticlesToSlices,
                                     SliceToHitCollection& outputSlicesToHits)
  {
    // Check for the special case in which there are no slices, and only the neutrino reconstruction was used on all hits
    if (settings.m_isNeutrinoRecoOnlyNoSlicing) {
      LArPandoraOutput::CopyAllHitsToSingleSlice(settings,
                                                 event,
                                                 instanceLabel,
                                                 pfoVector,
                                                 idToHitMap,
                                                 outputSlices,
                                                 outputParticlesToSlices,
                                                 outputSlicesToHits);
      return;
    }

    // Collect the slice pfos - one per slice (if there is no slicing instance, this vector will be empty)
    pandora::PfoVector slicePfos;
    LArPandoraOutput::GetPandoraSlices(pPrimaryPandora, slicePfos);

    // Make one slice per Pandora Slice pfo
    for (const pandora::ParticleFlowObject* const pSlicePfo : slicePfos)
      LArPandoraOutput::BuildSlice(
        pSlicePfo, event, instanceLabel, idToHitMap, outputSlices, outputSlicesToHits);

    // Make a slice for every remaining pfo hierarchy that wasn't already in a slice
    std::unordered_map<const pandora::ParticleFlowObject*, unsigned int> parentPfoToSliceIndexMap;
    for (unsigned int pfoId = 0; pfoId < pfoVector.size(); ++pfoId) {
      const pandora::ParticleFlowObject* const pPfo(pfoVector.at(pfoId));

      // If this PFO is the parent of a hierarchy we have yet to use, then add a new slice
      if (LArPandoraOutput::IsFromSlice(pPfo)) continue;

      if (lar_content::LArPfoHelper::GetParentPfo(pPfo) != pPfo) continue;

      if (!parentPfoToSliceIndexMap
             .emplace(pPfo,
                      LArPandoraOutput::BuildSlice(
                        pPfo, event, instanceLabel, idToHitMap, outputSlices, outputSlicesToHits))
             .second)
        throw cet::exception("LArPandora")
          << " LArPandoraOutput::BuildSlices --- found repeated primary particles ";
    }

    // Add the associations from PFOs to slices
    for (unsigned int pfoId = 0; pfoId < pfoVector.size(); ++pfoId) {
      const pandora::ParticleFlowObject* const pPfo(pfoVector.at(pfoId));

      // For PFOs that are from a Pandora slice, add the association and move on to the next PFO
      if (LArPandoraOutput::IsFromSlice(pPfo)) {
        LArPandoraOutput::AddAssociation(event,
                                         instanceLabel,
                                         pfoId,
                                         LArPandoraOutput::GetSliceIndex(pPfo),
                                         outputParticlesToSlices);
        continue;
      }

      // Get the parent of the particle
      const pandora::ParticleFlowObject* const pParent(
        lar_content::LArPfoHelper::GetParentPfo(pPfo));
      if (parentPfoToSliceIndexMap.find(pParent) == parentPfoToSliceIndexMap.end())
        throw cet::exception("LArPandora")
          << " LArPandoraOutput::BuildSlices --- found pfo without a parent in the input list ";

      // Add the association from the PFO to the slice
      LArPandoraOutput::AddAssociation(
        event, instanceLabel, pfoId, parentPfoToSliceIndexMap.at(pParent), outputParticlesToSlices);
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  unsigned int LArPandoraOutput::BuildDummySlice(SliceCollection& outputSlices)
  {
    // Make a slice with dummy properties
    const float bogusFloat(std::numeric_limits<float>::max());
    const recob::tracking::Point_t bogusPoint(bogusFloat, bogusFloat, bogusFloat);
    const recob::tracking::Vector_t bogusVector(bogusFloat, bogusFloat, bogusFloat);

    const unsigned int sliceIndex(outputSlices->size());
    outputSlices->emplace_back(
      sliceIndex, bogusPoint, bogusVector, bogusPoint, bogusPoint, bogusFloat, bogusFloat);

    return sliceIndex;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraOutput::CopyAllHitsToSingleSlice(
    const Settings& settings,
    const art::Event& event,
    const std::string& instanceLabel,
    const pandora::PfoVector& pfoVector,
    const IdToHitMap& idToHitMap,
    SliceCollection& outputSlices,
    PFParticleToSliceCollection& outputParticlesToSlices,
    SliceToHitCollection& outputSlicesToHits)
  {
    const unsigned int sliceIndex(LArPandoraOutput::BuildDummySlice(outputSlices));

    // Add all of the hits in the events to the slice
    HitVector hits;
    LArPandoraHelper::CollectHits(event, settings.m_hitfinderModuleLabel, hits);
    LArPandoraOutput::AddAssociation(event, instanceLabel, sliceIndex, hits, outputSlicesToHits);

    mf::LogDebug("LArPandora") << "Finding hits with label: " << settings.m_hitfinderModuleLabel
                               << std::endl;
    mf::LogDebug("LArPandora") << " - Found " << hits.size() << std::endl;
    mf::LogDebug("LArPandora") << " - Making associations " << outputSlicesToHits->size()
                               << std::endl;

    // Add all of the PFOs to the slice
    for (unsigned int pfoId = 0; pfoId < pfoVector.size(); ++pfoId)
      LArPandoraOutput::AddAssociation(
        event, instanceLabel, pfoId, sliceIndex, outputParticlesToSlices);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  unsigned int LArPandoraOutput::BuildSlice(const pandora::ParticleFlowObject* const pParentPfo,
                                            const art::Event& event,
                                            const std::string& instanceLabel,
                                            const IdToHitMap& idToHitMap,
                                            SliceCollection& outputSlices,
                                            SliceToHitCollection& outputSlicesToHits)
  {
    const unsigned int sliceIndex(LArPandoraOutput::BuildDummySlice(outputSlices));

    // Collect the pfos connected to the input primary pfos
    pandora::PfoList pfosInSlice;
    lar_content::LArPfoHelper::GetAllConnectedPfos(pParentPfo, pfosInSlice);
    pfosInSlice.sort(lar_content::LArPfoHelper::SortByNHits);

    // Collect the hits from the pfos in all views
    pandora::CaloHitList hits;
    for (const pandora::ParticleFlowObject* const pPfo : pfosInSlice) {
      for (const pandora::HitType& hitType :
           {pandora::TPC_VIEW_U, pandora::TPC_VIEW_V, pandora::TPC_VIEW_W}) {
        lar_content::LArPfoHelper::GetCaloHits(pPfo, hitType, hits);
        lar_content::LArPfoHelper::GetIsolatedCaloHits(pPfo, hitType, hits);
      }
    }

    // Add the associations to the hits
    for (const pandora::CaloHit* const pCaloHit : hits)
      LArPandoraOutput::AddAssociation(event,
                                       instanceLabel,
                                       sliceIndex,
                                       {LArPandoraOutput::GetHit(idToHitMap, pCaloHit)},
                                       outputSlicesToHits);

    return sliceIndex;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraOutput::BuildT0s(const art::Event& event,
                                  const std::string& instanceLabel,
                                  const pandora::PfoVector& pfoVector,
                                  T0Collection& outputT0s,
                                  PFParticleToT0Collection& outputParticlesToT0s)
  {
    size_t nextT0Id(0);
    for (unsigned int pfoId = 0; pfoId < pfoVector.size(); ++pfoId) {
      const pandora::ParticleFlowObject* const pPfo(pfoVector.at(pfoId));

      anab::T0 t0;
      if (!LArPandoraOutput::BuildT0(event, pPfo, pfoVector, nextT0Id, t0)) continue;

      LArPandoraOutput::AddAssociation(
        event, instanceLabel, pfoId, nextT0Id - 1, outputParticlesToT0s);
      outputT0s->push_back(t0);
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  recob::PFParticle LArPandoraOutput::BuildPFParticle(const pandora::ParticleFlowObject* const pPfo,
                                                      const size_t pfoId,
                                                      const pandora::PfoVector& pfoVector)
  {
    // Get parent Pfo ID
    const pandora::PfoList& parentList(pPfo->GetParentPfoList());
    if (parentList.size() > 1)
      throw cet::exception("LArPandora")
        << " LArPandoraOutput::BuildPFParticle --- this pfo has multiple parent particles ";

    const size_t parentId(parentList.empty() ?
                            recob::PFParticle::kPFParticlePrimary :
                            LArPandoraOutput::GetId(parentList.front(), pfoVector));

    // Get daughters Pfo IDs
    std::vector<size_t> daughterIds;
    for (const pandora::ParticleFlowObject* const pDaughterPfo : pPfo->GetDaughterPfoList())
      daughterIds.push_back(LArPandoraOutput::GetId(pDaughterPfo, pfoVector));

    std::sort(daughterIds.begin(), daughterIds.end());

    return recob::PFParticle(pPfo->GetParticleId(), pfoId, parentId, daughterIds);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  recob::Vertex LArPandoraOutput::BuildVertex(const pandora::Vertex* const pVertex,
                                              const size_t vertexId)
  {
    double pos[3] = {
      pVertex->GetPosition().GetX(), pVertex->GetPosition().GetY(), pVertex->GetPosition().GetZ()};
    return recob::Vertex(pos, vertexId);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraOutput::GetHitsInCluster(const pandora::Cluster* const pCluster,
                                          pandora::CaloHitVector& sortedHits)
  {
    if (!sortedHits.empty())
      throw cet::exception("LArPandora")
        << " LArPandoraOutput::GetHitsInCluster --- vector to hold hits is not empty ";

    pandora::CaloHitList hitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(hitList);
    hitList.insert(hitList.end(),
                   pCluster->GetIsolatedCaloHitList().begin(),
                   pCluster->GetIsolatedCaloHitList().end());

    sortedHits.insert(sortedHits.end(), hitList.begin(), hitList.end());
    std::sort(
      sortedHits.begin(), sortedHits.end(), lar_content::LArClusterHelper::SortHitsByPosition);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  std::vector<recob::Cluster> LArPandoraOutput::BuildClusters(
    util::GeometryUtilities const& gser,
    const pandora::Cluster* const pCluster,
    const pandora::ClusterList& clusterList,
    const CaloHitToArtHitMap& pandoraHitToArtHitMap,
    IdToIdVectorMap& pandoraClusterToArtClustersMap,
    std::vector<HitVector>& hitVectors,
    size_t& nextId,
    cluster::ClusterParamsAlgBase& algo)
  {
    std::vector<recob::Cluster> clusters;

    // Get the cluster ID and set up the map entry
    const size_t clusterId(LArPandoraOutput::GetId(pCluster, clusterList));
    if (!pandoraClusterToArtClustersMap.insert(IdToIdVectorMap::value_type(clusterId, {})).second)
      throw cet::exception("LArPandora")
        << " LArPandoraOutput::BuildClusters --- repeated clusters in input list ";

    pandora::CaloHitVector sortedHits;
    LArPandoraOutput::GetHitsInCluster(pCluster, sortedHits);

    HitArray hitArray; // hits organised by drift volume
    HitList isolatedHits;

    for (const pandora::CaloHit* const pCaloHit2D : sortedHits) {
      CaloHitToArtHitMap::const_iterator it(pandoraHitToArtHitMap.find(pCaloHit2D));
      if (it == pandoraHitToArtHitMap.end())
        throw cet::exception("LArPandora")
          << " LArPandoraOutput::BuildClusters --- couldn't find art hit for input pandora hit ";

      const art::Ptr<recob::Hit> hit(it->second);

      const geo::WireID wireID(hit->WireID());
      const unsigned int volID(100000 * wireID.Cryostat + wireID.TPC);
      hitArray[volID].push_back(hit);

      if (pCaloHit2D->IsIsolated()) isolatedHits.insert(hit);
    }

    if (hitArray.empty())
      throw cet::exception("LArPandora")
        << " LArPandoraOutput::BuildClusters --- found a cluster with no hits ";

    for (const HitArray::value_type& hitArrayEntry : hitArray) {
      const HitVector& clusterHits(hitArrayEntry.second);

      clusters.push_back(
        LArPandoraOutput::BuildCluster(gser, nextId, clusterHits, isolatedHits, algo));
      hitVectors.push_back(clusterHits);
      pandoraClusterToArtClustersMap.at(clusterId).push_back(nextId);

      nextId++;
    }

    return clusters;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  recob::Cluster LArPandoraOutput::BuildCluster(util::GeometryUtilities const& gser,
                                                const size_t id,
                                                const HitVector& hitVector,
                                                const HitList& isolatedHits,
                                                cluster::ClusterParamsAlgBase& algo)
  {
    if (hitVector.empty())
      throw cet::exception("LArPandora")
        << " LArPandoraOutput::BuildCluster --- No input hits were provided ";

    // Fill list of cluster properties
    geo::View_t view(geo::kUnknown);
    geo::PlaneID planeID;

    double startWire(+std::numeric_limits<float>::max()), sigmaStartWire(0.0);
    double startTime(+std::numeric_limits<float>::max()), sigmaStartTime(0.0);
    double endWire(-std::numeric_limits<float>::max()), sigmaEndWire(0.0);
    double endTime(-std::numeric_limits<float>::max()), sigmaEndTime(0.0);

    std::vector<recob::Hit const*> hits_for_params;
    hits_for_params.reserve(hitVector.size());

    for (const art::Ptr<recob::Hit>& hit : hitVector) {
      const double thisWire(hit->WireID().Wire);
      const double thisWireSigma(0.5);
      const double thisTime(hit->PeakTime());
      const double thisTimeSigma(double(2. * hit->RMS()));
      const geo::View_t thisView(hit->View());
      const geo::PlaneID thisPlaneID(hit->WireID().planeID());

      if (geo::kUnknown == view) {
        view = thisView;
        planeID = thisPlaneID;
      }

      if (!(thisView == view && thisPlaneID == planeID)) {
        throw cet::exception("LArPandora")
          << " LArPandoraOutput::BuildCluster --- Input hits have inconsistent plane IDs ";
      }

      hits_for_params.push_back(&*hit);

      if (isolatedHits.count(hit)) continue;

      if (thisWire < startWire || (thisWire == startWire && thisTime < startTime)) {
        startWire = thisWire;
        sigmaStartWire = thisWireSigma;
        startTime = thisTime;
        sigmaStartTime = thisTimeSigma;
      }

      if (thisWire > endWire || (thisWire == endWire && thisTime > endTime)) {
        endWire = thisWire;
        sigmaEndWire = thisWireSigma;
        endTime = thisTime;
        sigmaEndTime = thisTimeSigma;
      }
    }

    // feed the algorithm with all the cluster hits
    algo.SetHitsFromPointers(gser, hits_for_params);

    // create the recob::Cluster directly in the vector
    return cluster::ClusterCreator(gser,
                                   algo,                  // algo
                                   startWire,             // start_wire
                                   sigmaStartWire,        // sigma_start_wire
                                   startTime,             // start_tick
                                   sigmaStartTime,        // sigma_start_tick
                                   endWire,               // end_wire
                                   sigmaEndWire,          // sigma_end_wire
                                   endTime,               // end_tick
                                   sigmaEndTime,          // sigma_end_tick
                                   id,                    // ID
                                   view,                  // view
                                   planeID,               // plane
                                   recob::Cluster::Sentry // sentry
                                   )
      .move();
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  recob::SpacePoint LArPandoraOutput::BuildSpacePoint(const pandora::CaloHit* const pCaloHit,
                                                      const size_t spacePointId)
  {
    if (pandora::TPC_3D != pCaloHit->GetHitType())
      throw cet::exception("LArPandora")
        << " LArPandoraOutput::BuildSpacePoint --- trying to build a space point from a 2D hit";

    const pandora::CartesianVector point(pCaloHit->GetPositionVector());
    double xyz[3] = {point.GetX(), point.GetY(), point.GetZ()};

    // ATTN using dummy information
    double dxdydz[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // TODO: Fill in the error matrix
    double chi2(0.0);

    return recob::SpacePoint(xyz, dxdydz, chi2, spacePointId);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  bool LArPandoraOutput::BuildT0(const art::Event& e,
                                 const pandora::ParticleFlowObject* const pPfo,
                                 const pandora::PfoVector& pfoVector,
                                 size_t& nextId,
                                 anab::T0& t0)
  {
    const pandora::ParticleFlowObject* const pParent(lar_content::LArPfoHelper::GetParentPfo(pPfo));
    const float x0(pParent->GetPropertiesMap().count("X0") ? pParent->GetPropertiesMap().at("X0") :
                                                             0.f);

    auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
    auto const det_prop =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clock_data);
    const double cm_per_tick(det_prop.GetXTicksCoefficient());
    const double ns_per_tick(sampling_rate(clock_data));

    // ATTN: T0 values are currently calculated in nanoseconds relative to the trigger offset. Only non-zero values are outputted.
    const double T0(x0 * ns_per_tick / cm_per_tick);

    if (std::fabs(T0) <= std::numeric_limits<double>::epsilon()) return false;

    // Output T0 objects [arguments are:  time (nanoseconds);  trigger type (3 for TPC stitching!);  pfparticle SelfID code;  T0 ID code]
    t0 = anab::T0(T0, 3, LArPandoraOutput::GetId(pPfo, pfoVector), nextId++);

    return true;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------------------------------------------------

  LArPandoraOutput::Settings::Settings()
    : m_pPrimaryPandora(nullptr)
    , m_shouldRunStitching(false)
    , m_shouldProduceAllOutcomes(false)
    , m_shouldProduceTestBeamInteractionVertices(false)
    , m_isNeutrinoRecoOnlyNoSlicing(false)
  {}

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraOutput::Settings::Validate() const
  {
    if (!m_pPrimaryPandora)
      throw cet::exception("LArPandora")
        << " LArPandoraOutput::Settings::Validate --- primary Pandora instance does not exist ";

    if (!m_shouldProduceAllOutcomes) return;

    if (m_allOutcomesInstanceLabel.empty())
      throw cet::exception("LArPandora")
        << " LArPandoraOutput::Settings::Validate --- all outcomes instance label not set ";
  }

} // namespace lar_pandora
