/**
 *  @file  larpandora/LArPandoraInterface/LArPandoraHelper.cxx
 *
 *  @brief helper function for LArPandoraInterface producer module
 *
 */

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/MVAOutput.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Wire.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "Objects/ParticleFlowObject.h"
#include "Pandora/PandoraInternal.h"
#include "Pandora/PdgTable.h"

#include <iostream>
#include <limits>

namespace lar_pandora {

  void LArPandoraHelper::CollectWires(const art::Event& evt,
                                      const std::string& label,
                                      WireVector& wireVector)
  {
    art::Handle<std::vector<recob::Wire>> theWires;
    evt.getByLabel(label, theWires);

    if (!theWires.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find wires... " << std::endl;
      return;
    }
    else {
      mf::LogDebug("LArPandora") << "  Found: " << theWires->size() << " Wires " << std::endl;
    }

    for (unsigned int i = 0; i < theWires->size(); ++i) {
      const art::Ptr<recob::Wire> wire(theWires, i);
      wireVector.push_back(wire);
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::CollectHits(const art::Event& evt,
                                     const std::string& label,
                                     HitVector& hitVector)
  {
    art::Handle<std::vector<recob::Hit>> theHits;
    evt.getByLabel(label, theHits);

    if (!theHits.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find hits... " << std::endl;
      return;
    }
    else {
      mf::LogDebug("LArPandora") << "  Found: " << theHits->size() << " Hits " << std::endl;
    }

    for (unsigned int i = 0; i < theHits->size(); ++i) {
      const art::Ptr<recob::Hit> hit(theHits, i);
      hitVector.push_back(hit);
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::CollectNuGraphHitLabels(const art::Event& evt,
                                                 const std::string& label,
                                                 HitToPred& hitToPred)
  {
    float pFilter, pSemantic;

    art::Handle<std::vector<recob::Hit>> theHits;
    evt.getByLabel(label, theHits);    

    if (!theHits.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find hits... " << std::endl;
      return;
    }
    else {
      mf::LogDebug("LArPandora") << "  Found: " << theHits->size() << " Hits " << std::endl;
    }

    art::Handle<std::vector<anab::FeatureVector<1>>> filterHandle;
    evt.getByLabel(label, filterHandle);

    if (!filterHandle.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find filter label... " << std::endl;
      return;
    }

    art::Handle<std::vector<anab::FeatureVector<5>>> semanticHandle;
    evt.getByLabel(label, semanticHandle);

    if (!semanticHandle.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find semantic label... " << std::endl;
      return;
    }

    for (unsigned int i = 0; i < theHits->size(); ++i) {
      const art::Ptr<recob::Hit> hit(theHits, i);

      // filter 
      pFilter = filterHandle->at(i).at(0);

      // semantic 
      // encoding: <first_semantic_category><second_semantic_category>.<confidence>
      // pSemantic = hitToSemanticAssoc.at(i)->at(0); // placeholder :) this is just the score of the MIP category
      pSemantic = semanticHandle->at(i).at(0); ///< placeholder :) this is just the score of the MIP category
                                               ///< to do, actually perform the encoding -- this was just to test the infrastructure
      std::cout << pFilter << "\t" << pSemantic << std::endl;

      hitToPred[hit] = std::make_pair(pFilter, pSemantic);
    }

  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::CollectPFParticles(const art::Event& evt,
                                            const std::string& label,
                                            PFParticleVector& particleVector)
  {
    art::Handle<std::vector<recob::PFParticle>> theParticles;
    evt.getByLabel(label, theParticles);

    if (!theParticles.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find particles... " << std::endl;
      return;
    }
    else {
      mf::LogDebug("LArPandora") << "  Found: " << theParticles->size() << " PFParticles "
                                 << std::endl;
    }

    for (unsigned int i = 0; i < theParticles->size(); ++i) {
      const art::Ptr<recob::PFParticle> particle(theParticles, i);
      particleVector.push_back(particle);
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::CollectSpacePoints(const art::Event& evt,
                                            const std::string& label,
                                            SpacePointVector& spacePointVector,
                                            SpacePointsToHits& spacePointsToHits)
  {
    HitsToSpacePoints hitsToSpacePoints;
    return LArPandoraHelper::CollectSpacePoints(
      evt, label, spacePointVector, spacePointsToHits, hitsToSpacePoints);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::CollectSpacePoints(const art::Event& evt,
                                            const std::string& label,
                                            SpacePointVector& spacePointVector,
                                            SpacePointsToHits& spacePointsToHits,
                                            HitsToSpacePoints& hitsToSpacePoints)
  {
    art::Handle<std::vector<recob::SpacePoint>> theSpacePoints;
    evt.getByLabel(label, theSpacePoints);

    if (!theSpacePoints.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find spacepoints... " << std::endl;
      return;
    }
    else {
      mf::LogDebug("LArPandora") << "  Found: " << theSpacePoints->size() << " SpacePoints "
                                 << std::endl;
    }

    art::FindOneP<recob::Hit> theHitAssns(theSpacePoints, evt, label);
    for (unsigned int i = 0; i < theSpacePoints->size(); ++i) {
      const art::Ptr<recob::SpacePoint> spacepoint(theSpacePoints, i);
      spacePointVector.push_back(spacepoint);
      const art::Ptr<recob::Hit> hit = theHitAssns.at(i);
      spacePointsToHits[spacepoint] = hit;
      hitsToSpacePoints[hit] = spacepoint;
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::CollectClusters(const art::Event& evt,
                                         const std::string& label,
                                         ClusterVector& clusterVector,
                                         ClustersToHits& clustersToHits)
  {
    art::Handle<std::vector<recob::Cluster>> theClusters;
    evt.getByLabel(label, theClusters);

    if (!theClusters.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find clusters... " << std::endl;
      return;
    }
    else {
      mf::LogDebug("LArPandora") << "  Found: " << theClusters->size() << " Clusters " << std::endl;
    }

    art::FindManyP<recob::Hit> theHitAssns(theClusters, evt, label);
    for (unsigned int i = 0; i < theClusters->size(); ++i) {
      const art::Ptr<recob::Cluster> cluster(theClusters, i);
      clusterVector.push_back(cluster);

      const std::vector<art::Ptr<recob::Hit>> hits = theHitAssns.at(i);
      for (unsigned int j = 0; j < hits.size(); ++j) {
        const art::Ptr<recob::Hit> hit = hits.at(j);
        clustersToHits[cluster].push_back(hit);
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::CollectPFParticles(const art::Event& evt,
                                            const std::string& label,
                                            PFParticleVector& particleVector,
                                            PFParticlesToSpacePoints& particlesToSpacePoints)
  {
    art::Handle<std::vector<recob::PFParticle>> theParticles;
    evt.getByLabel(label, theParticles);

    if (!theParticles.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find particles... " << std::endl;
      return;
    }
    else {
      mf::LogDebug("LArPandora") << "  Found: " << theParticles->size() << " PFParticles "
                                 << std::endl;
    }

    art::FindManyP<recob::SpacePoint> theSpacePointAssns(theParticles, evt, label);
    for (unsigned int i = 0; i < theParticles->size(); ++i) {
      const art::Ptr<recob::PFParticle> particle(theParticles, i);
      particleVector.push_back(particle);

      const std::vector<art::Ptr<recob::SpacePoint>> spacepoints = theSpacePointAssns.at(i);
      for (unsigned int j = 0; j < spacepoints.size(); ++j) {
        const art::Ptr<recob::SpacePoint> spacepoint = spacepoints.at(j);
        particlesToSpacePoints[particle].push_back(spacepoint);
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::CollectPFParticles(const art::Event& evt,
                                            const std::string& label,
                                            PFParticleVector& particleVector,
                                            PFParticlesToClusters& particlesToClusters)
  {
    art::Handle<std::vector<recob::PFParticle>> theParticles;
    evt.getByLabel(label, theParticles);

    if (!theParticles.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find particles... " << std::endl;
      return;
    }
    else {
      mf::LogDebug("LArPandora") << "  Found: " << theParticles->size() << " PFParticles "
                                 << std::endl;
    }

    art::FindManyP<recob::Cluster> theClusterAssns(theParticles, evt, label);
    for (unsigned int i = 0; i < theParticles->size(); ++i) {
      const art::Ptr<recob::PFParticle> particle(theParticles, i);
      particleVector.push_back(particle);

      const std::vector<art::Ptr<recob::Cluster>> clusters = theClusterAssns.at(i);
      for (unsigned int j = 0; j < clusters.size(); ++j) {
        const art::Ptr<recob::Cluster> cluster = clusters.at(j);
        particlesToClusters[particle].push_back(cluster);
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::CollectPFParticleMetadata(const art::Event& evt,
                                                   const std::string& label,
                                                   PFParticleVector& particleVector,
                                                   PFParticlesToMetadata& particlesToMetadata)
  {
    art::Handle<std::vector<recob::PFParticle>> theParticles;
    evt.getByLabel(label, theParticles);

    if (!theParticles.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find particles... " << std::endl;
      return;
    }
    else {
      mf::LogDebug("LArPandora") << "  Found: " << theParticles->size() << " PFParticles "
                                 << std::endl;
    }

    art::FindManyP<larpandoraobj::PFParticleMetadata> theMetadataAssns(theParticles, evt, label);
    for (unsigned int i = 0; i < theParticles->size(); ++i) {
      const art::Ptr<recob::PFParticle> particle(theParticles, i);
      particleVector.push_back(particle);

      const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfParticleMetadataList =
        theMetadataAssns.at(i);
      for (unsigned int j = 0; j < pfParticleMetadataList.size(); ++j) {
        const art::Ptr<larpandoraobj::PFParticleMetadata> pfParticleMetadata =
          pfParticleMetadataList.at(j);
        particlesToMetadata[particle].push_back(pfParticleMetadata);
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::CollectShowers(const art::Event& evt,
                                        const std::string& label,
                                        ShowerVector& showerVector,
                                        PFParticlesToShowers& particlesToShowers)
  {
    art::Handle<std::vector<recob::Shower>> theShowers;
    evt.getByLabel(label, theShowers);

    if (!theShowers.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find showers... " << std::endl;
      return;
    }
    else {
      mf::LogDebug("LArPandora") << "  Found: " << theShowers->size() << " Showers " << std::endl;
    }

    art::FindManyP<recob::PFParticle> theShowerAssns(theShowers, evt, label);
    for (unsigned int i = 0; i < theShowers->size(); ++i) {
      const art::Ptr<recob::Shower> shower(theShowers, i);
      showerVector.push_back(shower);

      const std::vector<art::Ptr<recob::PFParticle>> particles = theShowerAssns.at(i);
      for (unsigned int j = 0; j < particles.size(); ++j) {
        const art::Ptr<recob::PFParticle> particle = particles.at(j);
        particlesToShowers[particle].push_back(shower);
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::CollectTracks(const art::Event& evt,
                                       const std::string& label,
                                       TrackVector& trackVector,
                                       PFParticlesToTracks& particlesToTracks)
  {
    art::Handle<std::vector<recob::Track>> theTracks;
    evt.getByLabel(label, theTracks);

    if (!theTracks.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find tracks... " << std::endl;
      return;
    }
    else {
      mf::LogDebug("LArPandora") << "  Found: " << theTracks->size() << " Tracks " << std::endl;
    }

    art::FindManyP<recob::PFParticle> theTrackAssns(theTracks, evt, label);
    for (unsigned int i = 0; i < theTracks->size(); ++i) {
      const art::Ptr<recob::Track> track(theTracks, i);
      trackVector.push_back(track);

      const std::vector<art::Ptr<recob::PFParticle>> particles = theTrackAssns.at(i);
      for (unsigned int j = 0; j < particles.size(); ++j) {
        const art::Ptr<recob::PFParticle> particle = particles.at(j);
        particlesToTracks[particle].push_back(track);
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::CollectTracks(const art::Event& evt,
                                       const std::string& label,
                                       TrackVector& trackVector,
                                       TracksToHits& tracksToHits)
  {
    art::Handle<std::vector<recob::Track>> theTracks;
    evt.getByLabel(label, theTracks);

    if (!theTracks.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find tracks... " << std::endl;
      return;
    }
    else {
      mf::LogDebug("LArPandora") << "  Found: " << theTracks->size() << " Tracks " << std::endl;
    }

    art::FindManyP<recob::Hit> theHitAssns(theTracks, evt, label);
    for (unsigned int i = 0; i < theTracks->size(); ++i) {
      const art::Ptr<recob::Track> track(theTracks, i);
      trackVector.push_back(track);

      const std::vector<art::Ptr<recob::Hit>> hits = theHitAssns.at(i);
      for (unsigned int j = 0; j < hits.size(); ++j) {
        const art::Ptr<recob::Hit> hit = hits.at(j);
        tracksToHits[track].push_back(hit);
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::CollectShowers(const art::Event& evt,
                                        const std::string& label,
                                        ShowerVector& showerVector,
                                        ShowersToHits& showersToHits)
  {
    art::Handle<std::vector<recob::Shower>> theShowers;
    evt.getByLabel(label, theShowers);

    if (!theShowers.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find showers... " << std::endl;
      return;
    }
    else {
      mf::LogDebug("LArPandora") << "  Found: " << theShowers->size() << " Showers " << std::endl;
    }

    art::FindManyP<recob::Hit> theHitAssns(theShowers, evt, label);
    for (unsigned int i = 0; i < theShowers->size(); ++i) {
      const art::Ptr<recob::Shower> shower(theShowers, i);
      showerVector.push_back(shower);

      const std::vector<art::Ptr<recob::Hit>> hits = theHitAssns.at(i);
      for (unsigned int j = 0; j < hits.size(); ++j) {
        const art::Ptr<recob::Hit> hit = hits.at(j);
        showersToHits[shower].push_back(hit);
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::CollectSeeds(const art::Event& evt,
                                      const std::string& label,
                                      SeedVector& seedVector,
                                      PFParticlesToSeeds& particlesToSeeds)
  {
    art::Handle<std::vector<recob::Seed>> theSeeds;
    evt.getByLabel(label, theSeeds);

    if (!theSeeds.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find seeds... " << std::endl;
      return;
    }
    else {
      mf::LogDebug("LArPandora") << "  Found: " << theSeeds->size() << " Seeds " << std::endl;
    }

    art::FindManyP<recob::PFParticle> theSeedAssns(theSeeds, evt, label);
    for (unsigned int i = 0; i < theSeeds->size(); ++i) {
      const art::Ptr<recob::Seed> seed(theSeeds, i);
      seedVector.push_back(seed);

      const std::vector<art::Ptr<recob::PFParticle>> particles = theSeedAssns.at(i);
      for (unsigned int j = 0; j < particles.size(); ++j) {
        const art::Ptr<recob::PFParticle> particle = particles.at(j);
        particlesToSeeds[particle].push_back(seed);
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::CollectSeeds(const art::Event& evt,
                                      const std::string& label,
                                      SeedVector& seedVector,
                                      SeedsToHits& seedsToHits)
  {
    art::Handle<std::vector<recob::Seed>> theSeeds;
    evt.getByLabel(label, theSeeds);

    if (!theSeeds.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find seeds... " << std::endl;
      return;
    }
    else {
      mf::LogDebug("LArPandora") << "  Found: " << theSeeds->size() << " Seeds " << std::endl;
    }

    art::FindOneP<recob::Hit> theHitAssns(theSeeds, evt, label);

    if (!theHitAssns.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find seed associations... " << std::endl;
      return;
    }

    for (unsigned int i = 0; i < theSeeds->size(); ++i) {
      const art::Ptr<recob::Seed> seed(theSeeds, i);
      seedVector.push_back(seed);
      const art::Ptr<recob::Hit> hit = theHitAssns.at(i);
      seedsToHits[seed] = hit;
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::CollectVertices(const art::Event& evt,
                                         const std::string& label,
                                         VertexVector& vertexVector,
                                         PFParticlesToVertices& particlesToVertices)
  {
    art::Handle<std::vector<recob::Vertex>> theVertices;
    evt.getByLabel(label, theVertices);

    if (!theVertices.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find vertices... " << std::endl;
      return;
    }
    else {
      mf::LogDebug("LArPandora") << "  Found: " << theVertices->size() << " Vertices " << std::endl;
    }

    art::FindManyP<recob::PFParticle> theVerticesAssns(theVertices, evt, label);
    for (unsigned int i = 0; i < theVertices->size(); ++i) {
      const art::Ptr<recob::Vertex> vertex(theVertices, i);
      vertexVector.push_back(vertex);

      const std::vector<art::Ptr<recob::PFParticle>> particles = theVerticesAssns.at(i);
      for (unsigned int j = 0; j < particles.size(); ++j) {
        const art::Ptr<recob::PFParticle> particle = particles.at(j);
        particlesToVertices[particle].push_back(vertex);
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::BuildPFParticleHitMaps(
    const PFParticleVector& particleVector,
    const PFParticlesToSpacePoints& particlesToSpacePoints,
    const SpacePointsToHits& spacePointsToHits,
    PFParticlesToHits& particlesToHits,
    HitsToPFParticles& hitsToParticles,
    const DaughterMode daughterMode)
  {
    // Build mapping from particle to particle ID for parent/daughter navigation
    PFParticleMap particleMap;

    for (PFParticleVector::const_iterator iter1 = particleVector.begin(),
                                          iterEnd1 = particleVector.end();
         iter1 != iterEnd1;
         ++iter1) {
      const art::Ptr<recob::PFParticle> particle = *iter1;
      particleMap[particle->Self()] = particle;
    }

    // Loop over hits and build mapping between reconstructed final-state particles and reconstructed hits
    for (PFParticlesToSpacePoints::const_iterator iter1 = particlesToSpacePoints.begin(),
                                                  iterEnd1 = particlesToSpacePoints.end();
         iter1 != iterEnd1;
         ++iter1) {
      const art::Ptr<recob::PFParticle> thisParticle = iter1->first;
      const art::Ptr<recob::PFParticle> particle(
        (kAddDaughters == daughterMode) ?
          LArPandoraHelper::GetFinalStatePFParticle(particleMap, thisParticle) :
          thisParticle);

      if ((kIgnoreDaughters == daughterMode) &&
          !LArPandoraHelper::IsFinalState(particleMap, particle))
        continue;

      const SpacePointVector& spacePointVector = iter1->second;

      for (SpacePointVector::const_iterator iter2 = spacePointVector.begin(),
                                            iterEnd2 = spacePointVector.end();
           iter2 != iterEnd2;
           ++iter2) {
        const art::Ptr<recob::SpacePoint> spacepoint = *iter2;

        SpacePointsToHits::const_iterator iter3 = spacePointsToHits.find(spacepoint);
        if (spacePointsToHits.end() == iter3)
          throw cet::exception("LArPandora") << " PandoraCollector::BuildPFParticleHitMaps --- "
                                                "Found a space point without an associated hit ";

        const art::Ptr<recob::Hit> hit = iter3->second;

        particlesToHits[particle].push_back(hit);
        hitsToParticles[hit] = particle;
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::BuildPFParticleHitMaps(const PFParticleVector& particleVector,
                                                const PFParticlesToClusters& particlesToClusters,
                                                const ClustersToHits& clustersToHits,
                                                PFParticlesToHits& particlesToHits,
                                                HitsToPFParticles& hitsToParticles,
                                                const DaughterMode daughterMode)
  {
    // Build mapping from particle to particle ID for parent/daughter navigation
    PFParticleMap particleMap;

    for (PFParticleVector::const_iterator iter1 = particleVector.begin(),
                                          iterEnd1 = particleVector.end();
         iter1 != iterEnd1;
         ++iter1) {
      const art::Ptr<recob::PFParticle> particle = *iter1;
      particleMap[particle->Self()] = particle;
    }

    // Loop over hits and build mapping between reconstructed final-state particles and reconstructed hits
    for (PFParticlesToClusters::const_iterator iter1 = particlesToClusters.begin(),
                                               iterEnd1 = particlesToClusters.end();
         iter1 != iterEnd1;
         ++iter1) {
      const art::Ptr<recob::PFParticle> thisParticle = iter1->first;
      const art::Ptr<recob::PFParticle> particle(
        (kAddDaughters == daughterMode) ?
          LArPandoraHelper::GetFinalStatePFParticle(particleMap, thisParticle) :
          thisParticle);

      if ((kIgnoreDaughters == daughterMode) &&
          !LArPandoraHelper::IsFinalState(particleMap, particle))
        continue;

      const ClusterVector& clusterVector = iter1->second;
      for (ClusterVector::const_iterator iter2 = clusterVector.begin(),
                                         iterEnd2 = clusterVector.end();
           iter2 != iterEnd2;
           ++iter2) {
        const art::Ptr<recob::Cluster> cluster = *iter2;

        ClustersToHits::const_iterator iter3 = clustersToHits.find(cluster);
        if (clustersToHits.end() == iter3)
          throw cet::exception("LArPandora") << " PandoraCollector::BuildPFParticleHitMaps --- "
                                                "Found a space point without an associated hit ";

        const HitVector& hitVector = iter3->second;
        for (HitVector::const_iterator iter4 = hitVector.begin(), iterEnd4 = hitVector.end();
             iter4 != iterEnd4;
             ++iter4) {
          const art::Ptr<recob::Hit> hit = *iter4;

          particlesToHits[particle].push_back(hit);
          hitsToParticles[hit] = particle;
        }
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::BuildPFParticleHitMaps(const art::Event& evt,
                                                const std::string& label,
                                                PFParticlesToHits& particlesToHits,
                                                HitsToPFParticles& hitsToParticles,
                                                const DaughterMode daughterMode,
                                                const bool useClusters)
  {
    return LArPandoraHelper::BuildPFParticleHitMaps(
      evt, label, label, particlesToHits, hitsToParticles, daughterMode, useClusters);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::BuildPFParticleHitMaps(const art::Event& evt,
                                                const std::string& label_pfpart,
                                                const std::string& label_middle,
                                                PFParticlesToHits& particlesToHits,
                                                HitsToPFParticles& hitsToParticles,
                                                const DaughterMode daughterMode,
                                                const bool useClusters)
  {
    // Use intermediate clusters
    if (useClusters) {
      PFParticleVector particleVector;
      PFParticlesToClusters particlesToClusters;

      ClusterVector clusterVector;
      ClustersToHits clustersToHits;

      LArPandoraHelper::CollectPFParticles(evt, label_pfpart, particleVector, particlesToClusters);
      LArPandoraHelper::CollectClusters(evt, label_middle, clusterVector, clustersToHits);

      LArPandoraHelper::BuildPFParticleHitMaps(particleVector,
                                               particlesToClusters,
                                               clustersToHits,
                                               particlesToHits,
                                               hitsToParticles,
                                               daughterMode);
    }

    // Use intermediate space points
    else {
      PFParticleVector particleVector;
      PFParticlesToSpacePoints particlesToSpacePoints;

      SpacePointVector spacePointVector;
      SpacePointsToHits spacePointsToHits;

      LArPandoraHelper::CollectPFParticles(
        evt, label_pfpart, particleVector, particlesToSpacePoints);
      LArPandoraHelper::CollectSpacePoints(evt, label_middle, spacePointVector, spacePointsToHits);

      LArPandoraHelper::BuildPFParticleHitMaps(particleVector,
                                               particlesToSpacePoints,
                                               spacePointsToHits,
                                               particlesToHits,
                                               hitsToParticles,
                                               daughterMode);
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::SelectNeutrinoPFParticles(const PFParticleVector& inputParticles,
                                                   PFParticleVector& outputParticles)
  {
    for (PFParticleVector::const_iterator iter = inputParticles.begin(),
                                          iterEnd = inputParticles.end();
         iter != iterEnd;
         ++iter) {
      const art::Ptr<recob::PFParticle> particle = *iter;

      if (LArPandoraHelper::IsNeutrino(particle)) outputParticles.push_back(particle);
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::SelectFinalStatePFParticles(const PFParticleVector& inputParticles,
                                                     PFParticleVector& outputParticles)
  {
    // Build mapping from particle to particle ID for parent/daughter navigation
    PFParticleMap particleMap;

    for (PFParticleVector::const_iterator iter = inputParticles.begin(),
                                          iterEnd = inputParticles.end();
         iter != iterEnd;
         ++iter) {
      const art::Ptr<recob::PFParticle> particle = *iter;
      particleMap[particle->Self()] = particle;
    }

    // Select final-state particles
    for (PFParticleVector::const_iterator iter = inputParticles.begin(),
                                          iterEnd = inputParticles.end();
         iter != iterEnd;
         ++iter) {
      const art::Ptr<recob::PFParticle> particle = *iter;

      if (LArPandoraHelper::IsFinalState(particleMap, particle))
        outputParticles.push_back(particle);
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::CollectCosmicTags(const art::Event& evt,
                                           const std::string& label,
                                           CosmicTagVector& cosmicTagVector,
                                           TracksToCosmicTags& tracksToCosmicTags)
  {
    art::Handle<std::vector<anab::CosmicTag>> theCosmicTags;
    evt.getByLabel(label, theCosmicTags);

    if (theCosmicTags.isValid()) {
      art::FindOneP<recob::Track> theCosmicAssns(
        theCosmicTags, evt, label); // We assume there is one tag per algorithm
      for (unsigned int i = 0; i < theCosmicTags->size(); ++i) {
        const art::Ptr<anab::CosmicTag> cosmicTag(theCosmicTags, i);
        const art::Ptr<recob::Track> track = theCosmicAssns.at(i);
        tracksToCosmicTags[track].push_back(
          cosmicTag); // We assume there could be multiple algorithms
        cosmicTagVector.push_back(cosmicTag);
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::CollectT0s(const art::Event& evt,
                                    const std::string& label,
                                    T0Vector& t0Vector,
                                    PFParticlesToT0s& particlesToT0s)
  {
    art::Handle<std::vector<anab::T0>> theT0s;
    evt.getByLabel(label, theT0s);

    if (theT0s.isValid()) {
      art::FindManyP<recob::PFParticle> theAssns(theT0s, evt, label);
      for (unsigned int i = 0; i < theT0s->size(); ++i) {
        const art::Ptr<anab::T0> theT0(theT0s, i);
        t0Vector.push_back(theT0);

        const std::vector<art::Ptr<recob::PFParticle>> particles = theAssns.at(i);
        for (unsigned int j = 0; j < particles.size(); ++j) {
          const art::Ptr<recob::PFParticle> theParticle = particles.at(j);
          particlesToT0s[theParticle].push_back(
            theT0); // We assume there could be multiple T0s per PFParticle
        }
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::CollectSimChannels(const art::Event& evt,
                                            const std::string& label,
                                            SimChannelVector& simChannelVector,
                                            bool& areSimChannelsValid)
  {
    art::Handle<std::vector<sim::SimChannel>> theSimChannels;
    evt.getByLabel(label, theSimChannels);

    if (!theSimChannels.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find sim channels... " << std::endl;
      areSimChannelsValid = false;
      return;
    }
    else {
      mf::LogDebug("LArPandora") << "  Found: " << theSimChannels->size() << " SimChannels "
                                 << std::endl;
      areSimChannelsValid = true;
    }

    for (unsigned int i = 0; i < theSimChannels->size(); ++i) {
      const art::Ptr<sim::SimChannel> channel(theSimChannels, i);
      simChannelVector.push_back(channel);
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::CollectMCParticles(const art::Event& evt,
                                            const std::string& label,
                                            MCParticleVector& particleVector)
  {
    art::Handle<RawMCParticleVector> theParticles;
    evt.getByLabel(label, theParticles);

    if (!theParticles.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find MC particles... " << std::endl;
      return;
    }
    else {
      mf::LogDebug("LArPandora") << "  Found: " << theParticles->size() << " MC particles "
                                 << std::endl;
    }

    for (unsigned int i = 0; i < theParticles->size(); ++i) {
      const art::Ptr<simb::MCParticle> particle(theParticles, i);
      particleVector.push_back(particle);
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::CollectGeneratorMCParticles(const art::Event& evt,
                                                     const std::string& label,
                                                     RawMCParticleVector& particleVector)
  {
    art::Handle<std::vector<simb::MCTruth>> mcTruthBlocks;
    evt.getByLabel(label, mcTruthBlocks);

    if (!mcTruthBlocks.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find MC truth blocks from generator... "
                                 << std::endl;
      return;
    }
    else {
      mf::LogDebug("LArPandora") << "  Found: " << mcTruthBlocks->size() << " MC truth blocks "
                                 << std::endl;
    }

    if (mcTruthBlocks->size() < 1)
      throw cet::exception("LArPandora") << " PandoraCollector::CollectGeneratorMCParticles --- "
                                            "Unexpected number of MC truth blocks ";

    for (size_t b = 0; b < mcTruthBlocks->size(); ++b) {
      const art::Ptr<simb::MCTruth> mcTruth(mcTruthBlocks, b);
      for (int i = 0; i < mcTruth->NParticles(); ++i) {
        particleVector.push_back(mcTruth->GetParticle(i));
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::CollectMCParticles(const art::Event& evt,
                                            const std::string& label,
                                            MCTruthToMCParticles& truthToParticles,
                                            MCParticlesToMCTruth& particlesToTruth)
  {
    art::Handle<RawMCParticleVector> theParticles;
    evt.getByLabel(label, theParticles);

    if (!theParticles.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find MC particles... " << std::endl;
      return;
    }
    else {
      mf::LogDebug("LArPandora") << "  Found: " << theParticles->size() << " MC particles "
                                 << std::endl;
    }

    art::FindOneP<simb::MCTruth> theTruthAssns(theParticles, evt, label);

    for (unsigned int i = 0, iEnd = theParticles->size(); i < iEnd; ++i) {
      const art::Ptr<simb::MCParticle> particle(theParticles, i);
      const art::Ptr<simb::MCTruth> truth(theTruthAssns.at(i));
      truthToParticles[truth].push_back(particle);
      particlesToTruth[particle] = truth;
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::BuildMCParticleHitMaps(const art::Event& evt,
                                                const HitVector& hitVector,
                                                const SimChannelVector& simChannelVector,
                                                HitsToTrackIDEs& hitsToTrackIDEs)
  {
    auto const clock_data =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

    SimChannelMap simChannelMap;

    for (SimChannelVector::const_iterator iter = simChannelVector.begin(),
                                          iterEnd = simChannelVector.end();
         iter != iterEnd;
         ++iter) {
      const art::Ptr<sim::SimChannel> simChannel = *iter;
      simChannelMap.insert(SimChannelMap::value_type(simChannel->Channel(), simChannel));
    }

    for (HitVector::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end();
         iter != iterEnd;
         ++iter) {
      const art::Ptr<recob::Hit> hit = *iter;

      SimChannelMap::const_iterator sIter = simChannelMap.find(hit->Channel());
      if (simChannelMap.end() == sIter) continue; // Hit has no truth information [continue]

      // ATTN: Need to convert TDCtick (integer) to TDC (unsigned integer) before passing to simChannel
      const raw::TDCtick_t start_tick(clock_data.TPCTick2TDC(hit->PeakTimeMinusRMS()));
      const raw::TDCtick_t end_tick(clock_data.TPCTick2TDC(hit->PeakTimePlusRMS()));
      const unsigned int start_tdc((start_tick < 0) ? 0 : start_tick);
      const unsigned int end_tdc(end_tick);

      if (start_tdc > end_tdc) continue; // Hit undershoots the readout window [continue]

      const art::Ptr<sim::SimChannel> simChannel = sIter->second;
      const TrackIDEVector trackCollection(simChannel->TrackIDEs(start_tdc, end_tdc));

      if (trackCollection.empty()) continue; // Hit has no truth information [continue]

      for (unsigned int iTrack = 0, iTrackEnd = trackCollection.size(); iTrack < iTrackEnd;
           ++iTrack) {
        const sim::TrackIDE trackIDE = trackCollection.at(iTrack);
        hitsToTrackIDEs[hit].push_back(trackIDE);
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::BuildMCParticleHitMaps(const HitsToTrackIDEs& hitsToTrackIDEs,
                                                const MCTruthToMCParticles& truthToParticles,
                                                MCParticlesToHits& particlesToHits,
                                                HitsToMCParticles& hitsToParticles,
                                                const DaughterMode daughterMode)
  {
    // Build mapping between particles and track IDs for parent/daughter navigation
    MCParticleMap particleMap;

    for (MCTruthToMCParticles::const_iterator iter1 = truthToParticles.begin(),
                                              iterEnd1 = truthToParticles.end();
         iter1 != iterEnd1;
         ++iter1) {
      const MCParticleVector& particleVector = iter1->second;
      for (MCParticleVector::const_iterator iter2 = particleVector.begin(),
                                            iterEnd2 = particleVector.end();
           iter2 != iterEnd2;
           ++iter2) {
        const art::Ptr<simb::MCParticle> particle = *iter2;
        particleMap[particle->TrackId()] = particle;
      }
    }

    // Loop over hits and build mapping between reconstructed hits and true particles
    for (HitsToTrackIDEs::const_iterator iter1 = hitsToTrackIDEs.begin(),
                                         iterEnd1 = hitsToTrackIDEs.end();
         iter1 != iterEnd1;
         ++iter1) {
      const art::Ptr<recob::Hit> hit = iter1->first;
      const TrackIDEVector& trackCollection = iter1->second;

      int bestTrackID(-1);
      float bestEnergyFrac(0.f);

      for (TrackIDEVector::const_iterator iter2 = trackCollection.begin(),
                                          iterEnd2 = trackCollection.end();
           iter2 != iterEnd2;
           ++iter2) {
        const sim::TrackIDE& trackIDE = *iter2;
        const int trackID(std::abs(trackIDE.trackID)); // TODO: Find out why std::abs is needed
        const float energyFrac(trackIDE.energyFrac);

        if (energyFrac > bestEnergyFrac) {
          bestEnergyFrac = energyFrac;
          bestTrackID = trackID;
        }
      }

      if (bestTrackID >= 0) {
        MCParticleMap::const_iterator iter3 = particleMap.find(bestTrackID);
        if (particleMap.end() == iter3)
          throw cet::exception("LArPandora") << " PandoraCollector::BuildMCParticleHitMaps --- "
                                                "Found a track ID without an MC Particle ";

        try {
          const art::Ptr<simb::MCParticle> thisParticle = iter3->second;
          const art::Ptr<simb::MCParticle> primaryParticle(
            LArPandoraHelper::GetFinalStateMCParticle(particleMap, thisParticle));
          const art::Ptr<simb::MCParticle> selectedParticle(
            (kAddDaughters == daughterMode) ? primaryParticle : thisParticle);

          if ((kIgnoreDaughters == daughterMode) && (selectedParticle != primaryParticle)) continue;

          if (!(LArPandoraHelper::IsVisible(selectedParticle))) continue;

          particlesToHits[selectedParticle].push_back(hit);
          hitsToParticles[hit] = selectedParticle;
        }
        catch (cet::exception& e) {
        }
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::BuildMCParticleHitMaps(const art::Event& evt,
                                                const std::string& label,
                                                const HitVector& hitVector,
                                                MCParticlesToHits& particlesToHits,
                                                HitsToMCParticles& hitsToParticles,
                                                const DaughterMode daughterMode)
  {
    SimChannelVector simChannelVector;
    MCTruthToMCParticles truthToParticles;
    MCParticlesToMCTruth particlesToTruth;
    HitsToTrackIDEs hitsToTrackIDEs;

    bool areSimChannelsValid(false);
    LArPandoraHelper::CollectSimChannels(evt, label, simChannelVector, areSimChannelsValid);

    LArPandoraHelper::CollectMCParticles(evt, label, truthToParticles, particlesToTruth);
    LArPandoraHelper::BuildMCParticleHitMaps(evt, hitVector, simChannelVector, hitsToTrackIDEs);
    LArPandoraHelper::BuildMCParticleHitMaps(
      hitsToTrackIDEs, truthToParticles, particlesToHits, hitsToParticles, daughterMode);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::BuildMCParticleHitMaps(const art::Event& evt,
                                                const std::string& hitLabel,
                                                const std::string& backtrackLabel,
                                                HitsToTrackIDEs& hitsToTrackIDEs)
  {
    // Start by getting the collection of Hits
    art::Handle<std::vector<recob::Hit>> theHits;
    evt.getByLabel(hitLabel, theHits);

    if (!theHits.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find hits... " << std::endl;
      return;
    }

    HitVector hitVector;

    for (unsigned int i = 0; i < theHits->size(); ++i) {
      const art::Ptr<recob::Hit> hit(theHits, i);
      hitVector.push_back(hit);
    }

    // Now get the associations between Hits and MCParticles
    std::vector<anab::BackTrackerHitMatchingData const*> backtrackerVector;

    MCParticleVector particleVector;

    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> particles_per_hit(
      theHits, evt, backtrackLabel);

    if (!particles_per_hit.isValid()) {
      mf::LogDebug("LArPandora") << "  Failed to find reco-truth matching... " << std::endl;
      return;
    }

    // Now loop over the hits and build a collection of IDEs
    for (HitVector::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end();
         iter != iterEnd;
         ++iter) {
      const art::Ptr<recob::Hit> hit = *iter;

      particleVector.clear();
      backtrackerVector.clear();
      particles_per_hit.get(hit.key(), particleVector, backtrackerVector);

      for (unsigned int j = 0; j < particleVector.size(); ++j) {
        const art::Ptr<simb::MCParticle> particle = particleVector[j];

        sim::TrackIDE trackIDE;
        trackIDE.trackID = particle->TrackId();
        trackIDE.energy = backtrackerVector[j]->energy;
        trackIDE.energyFrac = backtrackerVector[j]->ideFraction;

        hitsToTrackIDEs[hit].push_back(trackIDE);
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::BuildMCParticleHitMaps(const art::Event& evt,
                                                const std::string& truthLabel,
                                                const std::string& hitLabel,
                                                const std::string& backtrackLabel,
                                                MCParticlesToHits& particlesToHits,
                                                HitsToMCParticles& hitsToParticles,
                                                const DaughterMode daughterMode)
  {
    MCTruthToMCParticles truthToParticles;
    MCParticlesToMCTruth particlesToTruth;
    HitsToTrackIDEs hitsToTrackIDEs;

    LArPandoraHelper::CollectMCParticles(evt, truthLabel, truthToParticles, particlesToTruth);
    LArPandoraHelper::BuildMCParticleHitMaps(evt, hitLabel, backtrackLabel, hitsToTrackIDEs);
    LArPandoraHelper::BuildMCParticleHitMaps(
      hitsToTrackIDEs, truthToParticles, particlesToHits, hitsToParticles, daughterMode);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  template <typename T>
  void LArPandoraHelper::GetAssociatedHits(const art::Event& evt,
                                           const std::string& label,
                                           const std::vector<art::Ptr<T>>& inputVector,
                                           HitVector& associatedHits,
                                           const pandora::IntVector* const indexVector)
  {

    art::Handle<std::vector<T>> handle;
    evt.getByLabel(label, handle);
    art::FindManyP<recob::Hit> hitAssoc(handle, evt, label);

    if (indexVector != nullptr) {
      if (inputVector.size() != indexVector->size())
        throw cet::exception("LArPandora") << " PandoraHelper::GetAssociatedHits --- trying to use "
                                              "an index vector not matching input vector";

      // If indexVector is filled, sort hits according to trajectory points order
      for (int index : (*indexVector)) {
        const art::Ptr<T>& element = inputVector.at(index);
        const HitVector& hits = hitAssoc.at(element.key());
        associatedHits.insert(associatedHits.end(), hits.begin(), hits.end());
      }
    }
    else {
      // If indexVector is empty just loop through inputSpacePoints
      for (const art::Ptr<T>& element : inputVector) {
        const HitVector& hits = hitAssoc.at(element.key());
        associatedHits.insert(associatedHits.end(), hits.begin(), hits.end());
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::BuildMCParticleMap(const MCParticleVector& particleVector,
                                            MCParticleMap& particleMap)
  {
    for (MCParticleVector::const_iterator iter = particleVector.begin(),
                                          iterEnd = particleVector.end();
         iter != iterEnd;
         ++iter) {
      const art::Ptr<simb::MCParticle> particle = *iter;
      particleMap[particle->TrackId()] = particle;
      particleMap[particle->TrackId()] = particle;
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraHelper::BuildPFParticleMap(const PFParticleVector& particleVector,
                                            PFParticleMap& particleMap)
  {
    for (PFParticleVector::const_iterator iter = particleVector.begin(),
                                          iterEnd = particleVector.end();
         iter != iterEnd;
         ++iter) {
      const art::Ptr<recob::PFParticle> particle = *iter;
      particleMap[particle->Self()] = particle;
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  art::Ptr<recob::PFParticle> LArPandoraHelper::GetParentPFParticle(
    const PFParticleMap& particleMap,
    const art::Ptr<recob::PFParticle> inputParticle)
  {
    // Navigate upward through PFO daughter/parent links - return the top-level PF Particle
    int primaryTrackID(inputParticle->Self());

    if (!inputParticle->IsPrimary()) {
      while (1) {
        PFParticleMap::const_iterator pIter1 = particleMap.find(primaryTrackID);
        if (particleMap.end() == pIter1)
          throw cet::exception("LArPandora") << " PandoraCollector::GetParentPFParticle --- Found "
                                                "a PFParticle without a particle ID ";

        const art::Ptr<recob::PFParticle> primaryParticle = pIter1->second;
        if (primaryParticle->IsPrimary()) break;

        primaryTrackID = primaryParticle->Parent();
      }
    }

    PFParticleMap::const_iterator pIter2 = particleMap.find(primaryTrackID);
    if (particleMap.end() == pIter2)
      throw cet::exception("LArPandora")
        << " PandoraCollector::GetParentPFParticle --- Found a PFParticle without a particle ID ";

    const art::Ptr<recob::PFParticle> outputParticle = pIter2->second;
    return outputParticle;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  art::Ptr<recob::PFParticle> LArPandoraHelper::GetFinalStatePFParticle(
    const PFParticleMap& particleMap,
    const art::Ptr<recob::PFParticle> inputParticle)
  {
    // Navigate upward through PFO daughter/parent links - return the top-level non-neutrino PF Particle
    int primaryTrackID(inputParticle->Self());

    if (!inputParticle->IsPrimary()) {
      int parentTrackID(inputParticle->Parent());

      while (1) {
        PFParticleMap::const_iterator pIter1 = particleMap.find(parentTrackID);
        if (particleMap.end() == pIter1)
          throw cet::exception("LArPandora") << " PandoraCollector::GetFinalStatePFParticle --- "
                                                "Found a PFParticle without a particle ID ";

        const art::Ptr<recob::PFParticle> parentParticle = pIter1->second;
        if (LArPandoraHelper::IsNeutrino(parentParticle)) break;

        primaryTrackID = parentTrackID;

        if (parentParticle->IsPrimary()) break;

        parentTrackID = parentParticle->Parent();
      }
    }

    PFParticleMap::const_iterator pIter2 = particleMap.find(primaryTrackID);
    if (particleMap.end() == pIter2)
      throw cet::exception("LArPandora") << " PandoraCollector::GetFinalStatePFParticle --- Found "
                                            "a PFParticle without a particle ID ";

    const art::Ptr<recob::PFParticle> outputParticle = pIter2->second;
    return outputParticle;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  art::Ptr<simb::MCParticle> LArPandoraHelper::GetParentMCParticle(
    const MCParticleMap& particleMap,
    const art::Ptr<simb::MCParticle> inputParticle)
  {
    // Navigate upward through MC daughter/parent links - return the top-level MC particle
    int primaryTrackID(inputParticle->TrackId());
    int parentTrackID(inputParticle->Mother());

    while (1) {
      MCParticleMap::const_iterator pIter1 = particleMap.find(parentTrackID);
      if (particleMap.end() == pIter1) break; // Can't find MC Particle for this track ID [break]

      const art::Ptr<simb::MCParticle> particle = pIter1->second;

      primaryTrackID = parentTrackID;
      parentTrackID = particle->Mother();
    }

    MCParticleMap::const_iterator pIter2 = particleMap.find(primaryTrackID);
    if (particleMap.end() == pIter2)
      throw cet::exception("LArPandora")
        << " PandoraCollector::GetParentMCParticle --- Found a track ID without a MC particle ";

    const art::Ptr<simb::MCParticle> outputParticle = pIter2->second;
    return outputParticle;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  art::Ptr<simb::MCParticle> LArPandoraHelper::GetFinalStateMCParticle(
    const MCParticleMap& particleMap,
    const art::Ptr<simb::MCParticle> inputParticle)
  {
    // Navigate upward through MC daughter/parent links - collect this particle and all its parents
    MCParticleVector mcVector;

    int trackID(inputParticle->TrackId());

    while (1) {
      MCParticleMap::const_iterator pIter = particleMap.find(trackID);
      if (particleMap.end() == pIter) break; // Can't find MC Particle for this track ID [break]

      const art::Ptr<simb::MCParticle> particle = pIter->second;
      mcVector.push_back(particle);

      trackID = particle->Mother();
    }

    // Navigate downward through MC parent/daughter links - return the first long-lived charged particle
    for (MCParticleVector::const_reverse_iterator iter = mcVector.rbegin(),
                                                  iterEnd = mcVector.rend();
         iter != iterEnd;
         ++iter) {
      const art::Ptr<simb::MCParticle> nextParticle = *iter;

      if (LArPandoraHelper::IsVisible(nextParticle)) return nextParticle;
    }

    throw cet::exception("LArPandora"); // need to catch this exception
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  art::Ptr<recob::Track> LArPandoraHelper::GetPrimaryTrack(
    const PFParticlesToTracks& particlesToTracks,
    const art::Ptr<recob::PFParticle> particle)
  {
    PFParticlesToTracks::const_iterator tIter = particlesToTracks.find(particle);

    if (particlesToTracks.end() == tIter || tIter->second.empty())
      throw cet::exception("LArPandora")
        << " PandoraCollector::GetPrimaryTrack --- Failed to find associated track ";

    if (tIter->second.size() != 1)
      throw cet::exception("LArPandora")
        << " PandoraCollector::GetPrimaryTrack --- Found more than one associated track ";

    const art::Ptr<recob::Track> primaryTrack = *(tIter->second.begin());
    return primaryTrack;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  int LArPandoraHelper::GetGeneration(const PFParticleMap& particleMap,
                                      const art::Ptr<recob::PFParticle> inputParticle)
  {
    // Navigate upward through PFO daughter/parent links - return the top-level PF Particle
    int nGenerations(0);
    int primaryTrackID(inputParticle->Self());

    while (1) {
      PFParticleMap::const_iterator pIter = particleMap.find(primaryTrackID);
      if (particleMap.end() == pIter)
        throw cet::exception("LArPandora")
          << " PandoraCollector::GetGeneration --- Found a PFParticle without a particle ID ";

      ++nGenerations;

      const art::Ptr<recob::PFParticle> primaryParticle = pIter->second;
      if (primaryParticle->IsPrimary()) break;

      primaryTrackID = primaryParticle->Parent();
    }

    return nGenerations;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  int LArPandoraHelper::GetParentNeutrino(const PFParticleMap& particleMap,
                                          const art::Ptr<recob::PFParticle> daughterParticle)
  {
    art::Ptr<recob::PFParticle> parentParticle =
      LArPandoraHelper::GetParentPFParticle(particleMap, daughterParticle);

    if (LArPandoraHelper::IsNeutrino(parentParticle)) return parentParticle->PdgCode();

    if (parentParticle->IsPrimary()) return 0;

    const int parentID(parentParticle->Parent());

    PFParticleMap::const_iterator pIter = particleMap.find(parentID);
    if (particleMap.end() == pIter)
      throw cet::exception("LArPandora")
        << " PandoraCollector::GetParentNeutrino --- Found a PFParticle without a particle ID ";

    const art::Ptr<recob::PFParticle> neutrinoParticle = pIter->second;
    return neutrinoParticle->PdgCode();
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  bool LArPandoraHelper::IsFinalState(const PFParticleMap& particleMap,
                                      const art::Ptr<recob::PFParticle> daughterParticle)
  {
    if (LArPandoraHelper::IsNeutrino(daughterParticle)) return false;

    if (daughterParticle->IsPrimary()) return true;

    const int parentID(daughterParticle->Parent());

    PFParticleMap::const_iterator pIter = particleMap.find(parentID);
    if (particleMap.end() == pIter)
      throw cet::exception("LArPandora")
        << " PandoraCollector::IsFinalState --- Found a PFParticle without a particle ID ";

    const art::Ptr<recob::PFParticle> parentParticle = pIter->second;

    if (LArPandoraHelper::IsNeutrino(parentParticle)) return true;

    return false;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  bool LArPandoraHelper::IsNeutrino(const art::Ptr<recob::PFParticle> particle)
  {
    const int pdg(particle->PdgCode());

    // electron, muon, tau (use Pandora PDG tables)
    return ((pandora::NU_E == std::abs(pdg)) || (pandora::NU_MU == std::abs(pdg)) ||
            (pandora::NU_TAU == std::abs(pdg)));
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  bool LArPandoraHelper::IsTrack(const art::Ptr<recob::PFParticle> particle)
  {
    const int pdg(particle->PdgCode());

    // muon, pion, proton, kaon (use Pandora PDG tables)
    return ((pandora::MU_MINUS == std::abs(pdg)) || (pandora::PI_PLUS == std::abs(pdg)) ||
            (pandora::PROTON == std::abs(pdg)) || (pandora::K_PLUS == std::abs(pdg)));
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  bool LArPandoraHelper::IsShower(const art::Ptr<recob::PFParticle> particle)
  {
    const int pdg(particle->PdgCode());

    // electron, photon (use Pandora PDG tables)
    return ((pandora::E_MINUS == std::abs(pdg)) || (pandora::PHOTON == std::abs(pdg)));
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  bool LArPandoraHelper::IsVisible(const art::Ptr<simb::MCParticle> particle)
  {
    // Include long-lived charged particles
    const int pdg(particle->PdgCode());

    if ((pandora::E_MINUS == std::abs(pdg)) || (pandora::MU_MINUS == std::abs(pdg)) ||
        (pandora::PROTON == std::abs(pdg)) || (pandora::PI_PLUS == std::abs(pdg)) ||
        (pandora::K_PLUS == std::abs(pdg)) || (pandora::SIGMA_MINUS == std::abs(pdg)) ||
        (pandora::SIGMA_PLUS == std::abs(pdg)) || (pandora::HYPERON_MINUS == std::abs(pdg)) ||
        (pandora::PHOTON == std::abs(pdg)) || (pandora::NEUTRON == std::abs(pdg)))
      return true;

    // TODO: What about ions, neutrons, photons? (Have included neutrons and photons for now)

    return false;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  larpandoraobj::PFParticleMetadata LArPandoraHelper::GetPFParticleMetadata(
    const pandora::ParticleFlowObject* const pPfo)
  {
    return larpandoraobj::PFParticleMetadata(pPfo->GetPropertiesMap());
  }

  //------------------------------------------------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------------------------------------------------

  template void LArPandoraHelper::GetAssociatedHits(const art::Event&,
                                                    const std::string&,
                                                    const std::vector<art::Ptr<recob::Cluster>>&,
                                                    HitVector&,
                                                    const pandora::IntVector* const);

  template void LArPandoraHelper::GetAssociatedHits(const art::Event&,
                                                    const std::string&,
                                                    const std::vector<art::Ptr<recob::SpacePoint>>&,
                                                    HitVector&,
                                                    const pandora::IntVector* const);

} // namespace lar_pandora
