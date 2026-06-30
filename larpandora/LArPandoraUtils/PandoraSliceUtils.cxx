/**
*
* @file larpandora/LArPandoraUtils/PandoraSliceUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about Slices
*/

#include "larpandora/LArPandoraUtils/PandoraSliceUtils.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"

namespace lar_pandora
{

std::vector<art::Ptr<recob::Hit>> PandoraSliceUtils::GetHits(const art::Ptr<recob::Slice> &pSlice, const art::Event &evt, const std::string &label)
{    
    return PandoraSliceUtils::GetAssocProductVector<recob::Hit>(pSlice,evt,label,label);
}

} // namespace lar_pandora
