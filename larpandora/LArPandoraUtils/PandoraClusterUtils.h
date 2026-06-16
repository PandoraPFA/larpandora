/**
 *
 * @file dunereco/AnaUtils/PandoraClusterUtils.h
 *
 * @brief Utility containing helpful functions for end users to access information about Clusters
*/

#ifndef PANDORA_CLUSTER_UTILS_H
#define PANDORA_CLUSTER_UTILS_H

#include "art/Framework/Principal/Event.h"

#include "larpandora/LArPandoraUtils/PandoraUtilsBase.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"

#include <string>
#include <vector>

namespace lar_pandora
{
/**
 *
 * @brief PandoraClusterUtils class
 *
*/
class PandoraClusterUtils:PandoraUtilsBase
{
public:
    /**
    * @brief Get the hits associated with the cluster.
    *
    * @param cluster is the cluster for which we want the hits
    * @param evt is the underlying art event
    * @param label is the label for the cluster producer
    * 
    * @return vector of art::Ptrs to the hits 
    */
    static std::vector<art::Ptr<recob::Hit>> GetHits(const art::Ptr<recob::Cluster> &pCluster, const art::Event &evt, const std::string &label);
};

} // namespace dune_ana

#endif // PANDORA_CLUSTER_UTILS_H

