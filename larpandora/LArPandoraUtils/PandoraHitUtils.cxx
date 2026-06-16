/**
*
* @file dunereco/AnaUtils/PandoraHitUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about hits
*/

//STL
#include <algorithm>
//ROOT
//ART
#include "art/Framework/Services/Registry/ServiceHandle.h"
//LARSOFT
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
//DUNE
#include "larpandora/LArPandoraUtils/PandoraHitUtils.h"

namespace lar_pandora
{

std::vector<art::Ptr<recob::SpacePoint>> PandoraHitUtils::GetSpacePoints(const art::Ptr<recob::Hit> &pHit,
    const art::Event &evt, const std::string &hitLabel, const std::string &hitToSpacePointLabel)

{
    return PandoraHitUtils::GetAssocProductVector<recob::SpacePoint>(pHit,evt,hitLabel,hitToSpacePointLabel);
}

std::vector<art::Ptr<recob::Hit>> PandoraHitUtils::GetHitsOnPlane(const std::vector<art::Ptr<recob::Hit>> &hits, 
    const geo::PlaneID::PlaneID_t planeID)
{
    std::vector<art::Ptr<recob::Hit>> hitsOnPlane;
    std::vector<art::Ptr<recob::Hit>>::const_iterator hitIt = hits.begin();
    while ((hitIt 
        = std::find_if(hitIt, hits.end(), [planeID](const art::Ptr<recob::Hit> &hit){return hit->WireID().Plane==planeID; })) 
        != hits.end())
        hitsOnPlane.emplace_back(*(hitIt++));
    return hitsOnPlane;

}

double PandoraHitUtils::LifetimeCorrection(detinfo::DetectorClocksData const& clockData,
                                           detinfo::DetectorPropertiesData const& detProp,
                                           const art::Ptr<recob::Hit> &pHit)
{
    return PandoraHitUtils::LifetimeCorrection(clockData, detProp,
                                               pHit->PeakTime(), clockData.TriggerTime());
}

double PandoraHitUtils::LifetimeCorrection(detinfo::DetectorClocksData const& clockData,
                                           detinfo::DetectorPropertiesData const& detProp,
                                           const double timeInTicks, const double t0InMicroS)
{
    const double tpcSamplingRateInMicroS(clockData.TPCClock().TickPeriod());
    const double tpcTriggerOffsetInTicks(trigger_offset(clockData));

    const double timeCorrectedForT0((timeInTicks - tpcTriggerOffsetInTicks)*tpcSamplingRateInMicroS - t0InMicroS);
    const double tauLifetime(detProp.ElectronLifetime());

    return std::exp(timeCorrectedForT0/tauLifetime);
}

double PandoraHitUtils::LifetimeCorrectedTotalHitCharge(detinfo::DetectorClocksData const& clockData,
                                                        detinfo::DetectorPropertiesData const& detProp,
                                                        const std::vector<art::Ptr<recob::Hit> > &hits)
{
    double totalHitCharge(0);
    for (unsigned int iHit = 0; iHit < hits.size(); iHit++)
        totalHitCharge += hits[iHit]->Integral() * PandoraHitUtils::LifetimeCorrection(clockData, detProp, hits[iHit]);

    return totalHitCharge;
}

} // namespace lar_pandora
