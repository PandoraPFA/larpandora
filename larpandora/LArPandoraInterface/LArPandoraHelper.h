/**
 *  @file  larpandora/LArPandoraInterface/LArPandoraHelper.h
 *
 *  @brief helper function for LArPandoraInterface producer module
 *
 */
#ifndef LAR_PANDORA_HELPER_H
#define LAR_PANDORA_HELPER_H

#include "art/Framework/Principal/Event.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardataobj/Simulation/SimChannel.h"

#include <map>
#include <set>
#include <vector>

namespace anab {
  class CosmicTag;
  class T0;
}
namespace pandora {
  class ParticleFlowObject;
  class Vertex;
  typedef std::vector<int> IntVector;
}
namespace recob {
  class Cluster;
  class Hit;
  class PFParticle;
  class Seed;
  class Shower;
  class Slice;
  class SpacePoint;
  class Track;
  class Vertex;
  class Wire;
}
namespace larpandoraobj {
  class PFParticleMetadata;
}
namespace sim {
  class SimChannel;
  struct TrackIDE;
}
namespace simb {
  class MCParticle;
  class MCTruth;
}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora {

  typedef std::set<art::Ptr<recob::Hit>> HitList;

  typedef std::vector<art::Ptr<recob::Wire>> WireVector;
  typedef std::vector<art::Ptr<recob::Hit>> HitVector;
  typedef std::vector<art::Ptr<recob::SpacePoint>> SpacePointVector;
  typedef std::vector<art::Ptr<recob::Cluster>> ClusterVector;
  typedef std::vector<art::Ptr<recob::Seed>> SeedVector;
  typedef std::vector<art::Ptr<recob::Vertex>> VertexVector;
  typedef std::vector<art::Ptr<recob::Track>> TrackVector;
  typedef std::vector<art::Ptr<recob::Shower>> ShowerVector;
  typedef std::vector<art::Ptr<recob::PFParticle>> PFParticleVector;
  typedef std::vector<art::Ptr<simb::MCTruth>> MCTruthVector;
  typedef std::vector<art::Ptr<simb::MCParticle>> MCParticleVector;
  typedef std::vector<simb::MCParticle> RawMCParticleVector;
  typedef std::vector<art::Ptr<sim::SimChannel>> SimChannelVector;
  typedef std::vector<sim::TrackIDE> TrackIDEVector;
  typedef std::vector<art::Ptr<anab::CosmicTag>> CosmicTagVector;
  typedef std::vector<art::Ptr<anab::T0>> T0Vector;
  typedef std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> MetadataVector;

  typedef std::unordered_set<art::Ptr<recob::Hit>> HitSet;

  typedef std::map<art::Ptr<recob::PFParticle>, TrackVector> PFParticlesToTracks;
  typedef std::map<art::Ptr<recob::PFParticle>, ShowerVector> PFParticlesToShowers;
  typedef std::map<art::Ptr<recob::PFParticle>, ClusterVector> PFParticlesToClusters;
  typedef std::map<art::Ptr<recob::PFParticle>, SeedVector> PFParticlesToSeeds;
  typedef std::map<art::Ptr<recob::PFParticle>, VertexVector> PFParticlesToVertices;
  typedef std::map<art::Ptr<recob::PFParticle>, SpacePointVector> PFParticlesToSpacePoints;
  typedef std::map<art::Ptr<recob::PFParticle>, HitVector> PFParticlesToHits;
  typedef std::map<art::Ptr<recob::PFParticle>, MetadataVector> PFParticlesToMetadata;
  typedef std::map<art::Ptr<recob::Track>, HitVector> TracksToHits;
  typedef std::map<art::Ptr<recob::Shower>, HitVector> ShowersToHits;
  typedef std::map<art::Ptr<recob::Cluster>, HitVector> ClustersToHits;
  typedef std::map<art::Ptr<recob::Seed>, art::Ptr<recob::Hit>> SeedsToHits;
  typedef std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>> SpacePointsToHits;
  typedef std::map<art::Ptr<simb::MCTruth>, MCParticleVector> MCTruthToMCParticles;
  typedef std::map<art::Ptr<simb::MCTruth>, HitVector> MCTruthToHits;
  typedef std::map<art::Ptr<simb::MCTruth>, art::Ptr<recob::PFParticle>> MCTruthToPFParticles;
  typedef std::map<art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>> MCParticlesToMCTruth;
  typedef std::map<art::Ptr<simb::MCParticle>, HitVector> MCParticlesToHits;
  typedef std::map<art::Ptr<simb::MCParticle>, art::Ptr<recob::PFParticle>>
    MCParticlesToPFParticles;
  typedef std::map<art::Ptr<recob::Hit>, art::Ptr<recob::SpacePoint>> HitsToSpacePoints;
  typedef std::map<art::Ptr<recob::Hit>, art::Ptr<recob::PFParticle>> HitsToPFParticles;
  typedef std::map<art::Ptr<recob::Hit>, art::Ptr<simb::MCParticle>> HitsToMCParticles;
  typedef std::map<art::Ptr<recob::Hit>, art::Ptr<simb::MCTruth>> HitsToMCTruth;
  typedef std::map<art::Ptr<recob::Hit>, TrackIDEVector> HitsToTrackIDEs;
  typedef std::map<art::Ptr<recob::Track>, CosmicTagVector> TracksToCosmicTags;
  typedef std::map<art::Ptr<recob::PFParticle>, T0Vector> PFParticlesToT0s;

  typedef std::map<int, art::Ptr<recob::PFParticle>> PFParticleMap;
  typedef std::map<int, art::Ptr<recob::Cluster>> ClusterMap;
  typedef std::map<int, art::Ptr<recob::SpacePoint>> SpacePointMap;
  typedef std::map<int, art::Ptr<recob::Hit>> HitMap;
  typedef std::map<int, art::Ptr<simb::MCParticle>> MCParticleMap;
  typedef std::map<int, art::Ptr<sim::SimChannel>> SimChannelMap;

  typedef std::map<const pandora::ParticleFlowObject*, size_t> ThreeDParticleMap;
  typedef std::map<const pandora::Vertex*, unsigned int> ThreeDVertexMap;
  typedef std::map<int, HitVector> HitArray;

  /**
 *  @brief  LArPandoraHelper class
 */
  class LArPandoraHelper {
  public:
    /**
     *  @brief  DaughterMode enumeration
     */
    enum DaughterMode {
      kIgnoreDaughters = 0, // Only use parent particles
      kUseDaughters = 1,    // Use both parent and daughter partcles
      kAddDaughters = 2     // Absorb daughter particles into parent particles
    };

    /**
     *  @brief Collect the reconstructed wires from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the Wire list in the event
     *  @param wireVector the ouput vector of Wire objects
     */
    static void CollectWires(const art::Event& evt,
                             const std::string& label,
                             WireVector& wireVector);

    /**
     *  @brief Collect the reconstructed Hits from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the Hit list in the event
     *  @param hitVector the ouput vector of Hit objects
     */
    static void CollectHits(const art::Event& evt, const std::string& label, HitVector& hitVector);

    /**
     *  @brief Collect the reconstructed PFParticles from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the PFParticle list in the event
     *  @param particleVector the output vector of PFParticle objects
     */
    static void CollectPFParticles(const art::Event& evt,
                                   const std::string& label,
                                   PFParticleVector& particleVector);

    /**
     *  @brief Collect the reconstructed SpacePoints and associated hits from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the SpacePoint list in the event
     *  @param spacePointVector the output vector of SpacePoint objects
     *  @param spacePointsToHits the output map from SpacePoint to Hit objects
     */
    static void CollectSpacePoints(const art::Event& evt,
                                   const std::string& label,
                                   SpacePointVector& spacePointVector,
                                   SpacePointsToHits& spacePointsToHits);

    /**
     *  @brief Collect the reconstructed SpacePoints and associated hits from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the SpacePoint list in the event
     *  @param spacePointVector the output vector of SpacePoint objects
     *  @param spacePointsToHits the output map from SpacePoint to Hit objects
     *  @param hitsToSpacePoints the output map from Hit to SpacePoint objects
     */
    static void CollectSpacePoints(const art::Event& evt,
                                   const std::string& label,
                                   SpacePointVector& spacePointVector,
                                   SpacePointsToHits& spacePointsToHits,
                                   HitsToSpacePoints& hitsToSpacePoints);

    /**
     *  @brief Collect the reconstructed Clusters and associated hits from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the SpacePoint list in the event
     *  @param clusterVector the output vector of Cluster objects
     *  @param clustersToHits the output map from Cluster to Hit objects
     */
    static void CollectClusters(const art::Event& evt,
                                const std::string& label,
                                ClusterVector& clusterVector,
                                ClustersToHits& clustersToHits);

    /**
     *  @brief Collect the reconstructed PFParticles and associated SpacePoints from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the PFParticle list in the event
     *  @param particleVector the output vector of PFParticle objects
     *  @param particlesToSpacePoints the output map from PFParticle to SpacePoint objects
     */
    static void CollectPFParticles(const art::Event& evt,
                                   const std::string& label,
                                   PFParticleVector& particleVector,
                                   PFParticlesToSpacePoints& particlesToSpacePoints);

    /**
     *  @brief Collect the reconstructed PFParticles and associated Clusters from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the PFParticle list in the event
     *  @param particleVector the output vector of PFParticle objects
     *  @param particlesToClusters the output map from PFParticle to Cluster objects
     */
    static void CollectPFParticles(const art::Event& evt,
                                   const std::string& label,
                                   PFParticleVector& particleVector,
                                   PFParticlesToClusters& particlesToClusters);

    /**
     *  @brief Collect the reconstructed PFParticle Metadata from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the PFParticle list in the event
     *  @param particleVector the output vector of PFParticle objects
     *  @param particlesToSpacePoints the output map from PFParticle to PFParticleMetadata objects
     */
    static void CollectPFParticleMetadata(const art::Event& evt,
                                          const std::string& label,
                                          PFParticleVector& particleVector,
                                          PFParticlesToMetadata& particlesToMetadata);

    /**
     *  @brief Collect the reconstructed PFParticles and associated Showers from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the PFParticle list in the event
     *  @param showerVector the output vector of Shower objects
     *  @param particlesToShowers the output map from PFParticle to Shower objects
     */
    static void CollectShowers(const art::Event& evt,
                               const std::string& label,
                               ShowerVector& showerVector,
                               PFParticlesToShowers& particlesToShowers);

    /**
     *  @brief Collect the reconstructed Showers and associated Hits from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the PFParticle list in the event
     *  @param showerVector the output vector of Shower objects
     *  @param showersToHits the output map from Shower to Hit objects
     */
    static void CollectShowers(const art::Event& evt,
                               const std::string& label,
                               ShowerVector& showerVector,
                               ShowersToHits& showersToHits);

    /**
     *  @brief Collect the reconstructed PFParticles and associated Tracks from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the PFParticle list in the event
     *  @param trackVector the output vector of Track objects
     *  @param particlesToTracks the output map from PFParticle to Track objects
     */
    static void CollectTracks(const art::Event& evt,
                              const std::string& label,
                              TrackVector& trackVector,
                              PFParticlesToTracks& particlesToTracks);

    /**
     *  @brief Collect the reconstructed Tracks and associated Hits from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the PFParticle list in the event
     *  @param trackVector the output vector of Track objects
     *  @param tracksToHits the output map from Track to Hit objects
     */
    static void CollectTracks(const art::Event& evt,
                              const std::string& label,
                              TrackVector& trackVector,
                              TracksToHits& tracksToHits);

    /**
     *  @brief Collect the reconstructed PFParticles and associated Seeds from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the PFParticle list in the event
     *  @param seedVector the output vector of Seed objects
     *  @param particlesToSeeds the output map from PFParticle to Seed objects
     */
    static void CollectSeeds(const art::Event& evt,
                             const std::string& label,
                             SeedVector& seedVector,
                             PFParticlesToSeeds& particlesToSeeds);

    /**
     *  @brief Collect the reconstructed Seeds and associated Hits from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the PFParticle list in the event
     *  @param seedVector the output vector of Seed objects
     *  @param seedsToHits the output map from Seed to Hit objects
     */
    static void CollectSeeds(const art::Event& evt,
                             const std::string& label,
                             SeedVector& seedVector,
                             SeedsToHits& seedsToHits);

    /**
     *  @brief Collect the reconstructed PFParticles and associated Vertices from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the PFParticle list in the event
     *  @param vertexVector the output vector of Vertex objects
     *  @param particlesToVertices the output map from PFParticle to Vertex objects
     */
    static void CollectVertices(const art::Event& evt,
                                const std::string& label,
                                VertexVector& vertexVector,
                                PFParticlesToVertices& particlesToVertices);

    /**
     *  @brief Build mapping between PFParticles and Hits using PFParticle/SpacePoint/Hit maps
     *
     *  @param particleVector the input vector of PFParticle objects
     *  @param particlesToSpacePoints the input map from PFParticle to SpacePoint objects
     *  @param spacePointsToHits the input map from SpacePoint to Hit objects
     *  @param particlesToHits the output map from PFParticle to Hit objects
     *  @param hitsToParticles the output map from Hit to PFParticle objects
     *  @param daughterMode treatment of daughter particles in construction of maps
     */
    static void BuildPFParticleHitMaps(const PFParticleVector& particleVector,
                                       const PFParticlesToSpacePoints& particlesToSpacePoints,
                                       const SpacePointsToHits& spacePointsToHits,
                                       PFParticlesToHits& particlesToHits,
                                       HitsToPFParticles& hitsToParticles,
                                       const DaughterMode daughterMode = kUseDaughters);

    /**
     *  @brief Build mapping between PFParticles and Hits using PFParticle/Cluster/Hit maps
     *
     *  @param particleVector the input vector of PFParticle objects
     *  @param particlesToClusters the input map from PFParticle to Cluster objects
     *  @param clustersToHits the input map from Cluster to Hit objects
     *  @param particlesToHits the output map from PFParticle to Hit objects
     *  @param hitsToParticles the output map from Hit to PFParticle objects
     *  @param daughterMode treatment of daughter particles in construction of maps
     */
    static void BuildPFParticleHitMaps(const PFParticleVector& particleVector,
                                       const PFParticlesToClusters& particlesToClusters,
                                       const ClustersToHits& clustersToHits,
                                       PFParticlesToHits& particlesToHits,
                                       HitsToPFParticles& hitsToParticles,
                                       const DaughterMode daughterMode = kUseDaughters);

    /**
     *  @brief Build mapping between PFParticles and Hits starting from ART event record
     *
     *  @param evt the ART event record
     *  @param label_pfpart the label for the PFParticle list in the event
     *  @param label_space the label for the Intermediate list in the event
     *  @param particlesToHits output map from PFParticle to Hit objects
     *  @param hitsToParticles output map from Hit to PFParticle objects
     *  @param daughterMode treatment of daughter particles in construction of maps
     *  @param useClusters choice of intermediate object (true for Clusters, false for SpacePoints)
     */
    static void BuildPFParticleHitMaps(const art::Event& evt,
                                       const std::string& label_pfpart,
                                       const std::string& label_mid,
                                       PFParticlesToHits& particlesToHits,
                                       HitsToPFParticles& hitsToParticles,
                                       const DaughterMode daughterMode = kUseDaughters,
                                       const bool useClusters = true);

    /**
     *  @brief Build mapping between PFParticles and Hits starting from ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the PFParticle list in the event
     *  @param particlesToHits output map from PFParticle to Hit objects
     *  @param hitsToParticles output map from Hit to PFParticle objects
     *  @param daughterMode treatment of daughter particles in construction of maps
     *  @param useClusters choice of intermediate object (true for Clusters, false for SpacePoints)
     */
    static void BuildPFParticleHitMaps(const art::Event& evt,
                                       const std::string& label,
                                       PFParticlesToHits& particlesToHits,
                                       HitsToPFParticles& hitsToParticles,
                                       const DaughterMode daughterMode = kUseDaughters,
                                       const bool useClusters = true);

    /**
     *  @brief Collect a vector of cosmic tags from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the cosmic tag information in the event
     *  @param cosmicTagVector output vector of CosmicTag objects
     *  @param tracksToCosmicTags output map from tracks to cosmic tags
     */
    static void CollectCosmicTags(const art::Event& evt,
                                  const std::string& label,
                                  CosmicTagVector& cosmicTagVector,
                                  TracksToCosmicTags& tracksToCosmicTags);

    /**
     *  @brief Collect a vector of T0s from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the T0 information in the event
     *  @param t0Vector output vector of T0 objects
     *  @param particlesToT0s output map from PParticles to T0s
     */
    static void CollectT0s(const art::Event& evt,
                           const std::string& label,
                           T0Vector& t0Vector,
                           PFParticlesToT0s& particlesToT0s);

    /**
     *  @brief Collect a vector of SimChannel objects from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the truth information in the event
     *  @param simChannelVector output vector of SimChannel objects
     *  @param areSimChannelsValid boolean parameter indicating if the sim channel collection exists
     */
    static void CollectSimChannels(const art::Event& evt,
                                   const std::string& label,
                                   SimChannelVector& simChannelVector,
                                   bool& areSimChannelsValid);

    /**
     *  @brief Collect a vector of MCParticle objects from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the truth information in the event
     *  @param particleVector the output vector of MCParticle objects
     */
    static void CollectMCParticles(const art::Event& evt,
                                   const std::string& label,
                                   MCParticleVector& particleVector);

    /**
     *  @brief Collect a vector of MCParticle objects from the generator in the ART event record.  ATTN: This function is
     *         needed as accessing generator (opposed to Geant4) level MCParticles requires use of MCTruth block.
     *
     *  @param evt the ART event record
     *  @param label the label for the truth information in the generator
     *  @param particleVector the output vector of MCParticle objects
     */
    static void CollectGeneratorMCParticles(const art::Event& evt,
                                            const std::string& label,
                                            RawMCParticleVector& particleVector);

    /**
     *  @brief Collect truth information from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the truth information in the event
     *  @param truthToParticles output map from MCTruth to MCParticle objects
     *  @param particlesToTruth output map from MCParticle to MCTruth objects
     */
    static void CollectMCParticles(const art::Event& evt,
                                   const std::string& label,
                                   MCTruthToMCParticles& truthToParticles,
                                   MCParticlesToMCTruth& particlesToTruth);

    /**
     *  @brief Collect the links from reconstructed hits to their true energy deposits
     *
     *  @param evt the art event containers
     *  @param hitVector the input vector of reconstructed hits
     *  @param simChannelVector the input vector of SimChannels
     *  @param hitsToTrackIDEs the out map from hits to true energy deposits
     */
    static void BuildMCParticleHitMaps(const art::Event& evt,
                                       const HitVector& hitVector,
                                       const SimChannelVector& simChannelVector,
                                       HitsToTrackIDEs& hitsToTrackIDEs);

    /**
     *  @brief Build mapping between Hits and MCParticles, starting from Hit/TrackIDE/MCParticle information
     *
     *  @param hitsToTrackIDEs the input map from hits to true energy deposits
     *  @param truthToParticles the input map of truth information
     *  @param particlesToHits the mapping between true particles and reconstructed hits
     *  @param hitsToParticles the mapping between reconstructed hits and true particles
     *  @param daughterMode treatment of daughter particles in construction of maps
     */
    static void BuildMCParticleHitMaps(const HitsToTrackIDEs& hitsToTrackIDEs,
                                       const MCTruthToMCParticles& truthToParticles,
                                       MCParticlesToHits& particlesToHits,
                                       HitsToMCParticles& hitsToParticles,
                                       const DaughterMode daughterMode = kUseDaughters);

    /**
     *  @brief Build mapping between Hits and MCParticles, starting from ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the truth information in the event
     *  @param hitVector the input vector of reconstructed hits
     *  @param particlesToHits the output mapping between true particles and reconstructed hits
     *  @param hitsToParticles the output mapping between reconstructed hits and true particles
     *  @param daughterMode treatment of daughter particles in construction of maps
     */
    static void BuildMCParticleHitMaps(const art::Event& evt,
                                       const std::string& label,
                                       const HitVector& hitVector,
                                       MCParticlesToHits& particlesToHits,
                                       HitsToMCParticles& hitsToParticles,
                                       const DaughterMode daughterMode = kUseDaughters);

    /**
     *  @brief  Get mapping between hits and true energy deposits using back-tracker information
     *
     *  @param  evt the event record
     *  @param  hitLabel the label of the collection of hits
     *  @param  backtrackLabel the label of the collection of back-tracker information
     *  @param  hitsToTrackIDEs the output map between hits and true energy deposits
     */
    static void BuildMCParticleHitMaps(const art::Event& evt,
                                       const std::string& hitLabel,
                                       const std::string& backtrackLabel,
                                       HitsToTrackIDEs& hitsToTrackIDEs);

    /**
     *  @brief Build mapping between Hits and MCParticles, starting from Hit/TrackIDE/MCParticle information
     *
     *  @param evt the event record
     *  @param truthLabel the label describing the G4 truth information
     *  @param hitLabel the label describing the hit collection
     *  @param backtrackLabel the label describing the back-tracker information
     *  @param particlesToHits the mapping between true particles and reconstructed hits
     *  @param hitsToParticles the mapping between reconstructed hits and true particles
     *  @param daughterMode treatment of daughter particles in construction of maps
     */
    static void BuildMCParticleHitMaps(const art::Event& evt,
                                       const std::string& truthLabel,
                                       const std::string& hitLabel,
                                       const std::string& backtrackLabel,
                                       MCParticlesToHits& particlesToHits,
                                       HitsToMCParticles& hitsToParticles,
                                       const DaughterMode daughterMode = kUseDaughters);

    /**
     *  @brief  Get all hits associated with input clusters
     *
     *  @param  evt the event containing the hits
     *  @param  label the label of the collection producing PFParticles
     *  @param  input vector input of T (clusters, spacepoints)
     *  @param  associatedHits output hits associated with T
     *  @param  indexVector vector of spacepoint indices reflecting trajectory points sorting order
     */
    template <typename T>
    static void GetAssociatedHits(const art::Event& evt,
                                  const std::string& label,
                                  const std::vector<art::Ptr<T>>& inputVector,
                                  HitVector& associatedHits,
                                  const pandora::IntVector* const indexVector = nullptr);

    /**
     *  @brief Select reconstructed neutrino particles from a list of all reconstructed particles
     *
     *  @param inputParticles the input vector of all particles (it has to be all of them!)
     *  @param outputParticles the output vector of final-state particles
     */
    static void SelectNeutrinoPFParticles(const PFParticleVector& inputParticles,
                                          PFParticleVector& outputParticles);

    /**
     *  @brief Select final-state reconstructed particles from a list of all reconstructed particles
     *
     *  @param inputParticles the input vector of all particles (it has to be all of them!)
     *  @param outputParticles the output vector of final-state particles
     */
    static void SelectFinalStatePFParticles(const PFParticleVector& inputParticles,
                                            PFParticleVector& outputParticles);

    /**
     *  @brief Build particle maps for true particles
     *
     *  @param particleVector the input vector of true particles
     *  @param particleMap the output mapping between true particle and true track ID
     */
    static void BuildMCParticleMap(const MCParticleVector& particleVector,
                                   MCParticleMap& particleMap);

    /**
     *  @brief Build particle maps for reconstructed particles
     *
     *  @param particleVector the input vector of reconstructed particles
     *  @param particleMap the output mapping between reconstructed particles and particle ID
     */
    static void BuildPFParticleMap(const PFParticleVector& particleVector,
                                   PFParticleMap& particleMap);

    /**
     *  @brief Return the top-level parent particle by navigating up the chain of parent/daughter associations
     *
     *  @param particleMap the mapping between reconstructed particle and particle ID
     *  @param daughterParticle the input PF particle
     *
     *  @return the top-level parent particle
     */
    static art::Ptr<recob::PFParticle> GetParentPFParticle(
      const PFParticleMap& particleMap,
      const art::Ptr<recob::PFParticle> daughterParticle);

    /**
     *  @brief Return the final-state parent particle by navigating up the chain of parent/daughter associations
     *
     *  @param particleMap the mapping between reconstructed particle and particle ID
     *  @param daughterParticle the input PF particle
     *
     *  @return the final-state parent particle
     */
    static art::Ptr<recob::PFParticle> GetFinalStatePFParticle(
      const PFParticleMap& particleMap,
      const art::Ptr<recob::PFParticle> daughterParticle);

    /**
     *  @brief Return the top-level parent particle by navigating up the chain of parent/daughter associations
     *
     *  @param particleMap the mapping between true particle and true track ID
     *  @param daughterParticle the input MC particle
     *
     *  @return the top-level parent particle
     */
    static art::Ptr<simb::MCParticle> GetParentMCParticle(
      const MCParticleMap& particleMap,
      const art::Ptr<simb::MCParticle> daughterParticle);

    /**
     *  @brief Return the final-state parent particle by navigating up the chain of parent/daughter associations
     *
     *  @param particleMap the mapping between true particle and true track ID
     *  @param daughterParticle the input MC particle
     *
     *  @return the final-state parent particle
     */
    static art::Ptr<simb::MCParticle> GetFinalStateMCParticle(
      const MCParticleMap& particleMap,
      const art::Ptr<simb::MCParticle> daughterParticle);

    /**
     *  @brief Return the primary track associated with a PFParticle
     *
     *  @param particlesToTracks the mapping between particles and tracks
     *  @param particle  the input particle
     */
    static art::Ptr<recob::Track> GetPrimaryTrack(const PFParticlesToTracks& particlesToTracks,
                                                  const art::Ptr<recob::PFParticle> particle);

    /**
     *  @brief Return the generation of this particle (first generation if primary)
     *
     *  @param particleMap the mapping between reconstructed particle and particle ID
     *  @param daughterParticle the input daughter particle
     *
     *  @return the nth generation in the particle hierarchy
     */
    static int GetGeneration(const PFParticleMap& particleMap,
                             const art::Ptr<recob::PFParticle> daughterParticle);

    /**
     *  @brief Return the parent neutrino PDG code (or zero for cosmics) for a given reconstructed particle
     *
     *  @param particleMap the mapping between reconstructed particle and particle ID
     *  @param daughterParticle the input daughter particle
     *
     *  @return the PDG code of the parent neutrinos (or zero for cosmics)
     */
    static int GetParentNeutrino(const PFParticleMap& particleMap,
                                 const art::Ptr<recob::PFParticle> daughterParticle);

    /**
     *  @brief Determine whether a particle has been reconstructed as a final-state particle
     *
     *  @param particleMap the mapping between reconstructed particle and particle ID
     *  @param daughterParticle the input daughter particle
     *
     *  @return true/false
     */
    static bool IsFinalState(const PFParticleMap& particleMap,
                             const art::Ptr<recob::PFParticle> daughterParticle);

    /**
     *  @brief Determine whether a particle has been reconstructed as a neutrino
     *
     *  @param particle the input particle
     *
     *  @return true/false
     */
    static bool IsNeutrino(const art::Ptr<recob::PFParticle> particle);

    /**
     *  @brief Determine whether a particle has been reconstructed as track-like
     *
     *  @param particle the input particle
     *
     *  @return true/false
     */
    static bool IsTrack(const art::Ptr<recob::PFParticle> particle);

    /**
     *  @brief Determine whether a particle has been reconstructed as shower-like
     *
     *  @param particle the input particle
     *
     *  @return true/false
     */
    static bool IsShower(const art::Ptr<recob::PFParticle> particle);

    /**
     *  @brief  Determine whether a particle is visible (i.e. long-lived charged particle)
     *
     *  @param  particle the input mc particle
     *
     *  @return true/false
     */
    static bool IsVisible(const art::Ptr<simb::MCParticle> particle);

    /**
     *  @brief  Get metadata associated to a PFO
     *
     *  @param  pPfo input ParticleFlowObject
     *
     *  @return larpandoraobj::PFParticleMetadata
     */
    static larpandoraobj::PFParticleMetadata GetPFParticleMetadata(
      const pandora::ParticleFlowObject* const pPfo);

    /**
     *  @brief  Check if geometry has collection views only
     *
     *  @param  
     *
     *  @return true/false
     */
    static bool IsDualPhase();


  };

} // namespace lar_pandora

#endif //  LAR_PANDORA_HELPER_H
