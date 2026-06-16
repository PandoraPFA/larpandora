/**
*
* @file larpandora/LArPandoraUtils/PandoraSpacePointUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about SpacePoints
*/

#include "larpandora/LArPandoraUtils/PandoraSpacePointUtils.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"

namespace lar_pandora
{

std::vector<art::Ptr<recob::Hit>> PandoraSpacePointUtils::GetHits(const art::Ptr<recob::SpacePoint> &pSpacepoint, const art::Event &evt, const std::string &label)
{    
    return PandoraSpacePointUtils::GetAssocProductVector<recob::Hit>(pSpacepoint,evt,label,label);
}

} // namespace lar_pandora


