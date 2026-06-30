/**
 *
 * @file larpandora/LArPandoraUtils/PandoraSliceUtils.h
 *
 * @brief Utility containing helpful functions for end users to access information about Slices
*/

#ifndef PANDORA_SLICE_UTILS_H
#define PANDORA_SLICE_UTILS_H

#include "art/Framework/Principal/Event.h"

#include "larpandora/LArPandoraUtils/PandoraUtilsBase.h"

#include <string>
#include <vector>

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"

namespace lar_pandora
{
/**
 *
 * @brief DUNEAnaSliceUtils class
 *
*/
class PandoraSliceUtils:PandoraUtilsBase
{
public:
    /**
    * @brief Get the hits associated with the slice.
    *
    * @param slice is the slice for which we want the hits
    * @param evt is the underlying art event
    * @param label is the label for the slice producer
    * 
    * @return vector of art::Ptrs to the hits
    */
    static std::vector<art::Ptr<recob::Hit>> GetHits(const art::Ptr<recob::Slice> &pSlice, const art::Event &evt, const std::string &label);
};

} // namespace lar_pandora

#endif // PANDORA_SLICE_UTILS_H

