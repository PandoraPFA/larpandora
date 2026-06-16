/**
*
* @file dunereco/AnaUtils/PandoraClusterUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about Clusters
*/

#include "larpandora/LArPandoraUtils/PandoraClusterUtils.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"

namespace lar_pandora
{

std::vector<art::Ptr<recob::Hit>> PandoraClusterUtils::GetHits(const art::Ptr<recob::Cluster> &pCluster, const art::Event &evt, const std::string &label)
{    
    return PandoraClusterUtils::GetAssocProductVector<recob::Hit>(pCluster,evt,label,label);
}


} // namespace dune_ana


