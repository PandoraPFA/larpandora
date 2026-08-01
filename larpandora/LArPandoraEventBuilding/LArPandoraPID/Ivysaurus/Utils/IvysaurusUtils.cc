/**
 *  @file   larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Utils/IvysaurusUtils.cc
 *
 *  @brief  Helpers used by the Ivysaurus PID
 */
//LArSoft
#include "larcore/Geometry/WireReadout.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larpandora/LArPandoraUtils/PandoraHitUtils.h"
#include "larpandora/LArPandoraUtils/PandoraSpacePointUtils.h"
#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"
#include "larpandora/LArPandoraInterface/Detectors/GetDetectorType.h"
#include "larpandora/LArPandoraInterface/Detectors/LArPandoraDetectorType.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Utils/IvysaurusUtils.h"
#include "larsim/Utils/TruthMatchUtils.h"
//ROOT
#include "TVector3.h"

namespace IvysaurusUtils
{

void ObtainPandoraHitPositionAndWidth(const art::Event &evt, const art::Ptr<recob::Hit> hit, 
    const PandoraView hitType, TVector3 &pandoraHitPosition, float &width)
{
    auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout const>()->Get();
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

    const geo::WireID hitWireID = hit->WireID();
    const geo::CryostatID cryostatID(hitWireID.Cryostat);
    const double hitTime = hit->PeakTime();
    const double hitTimeStart = hit->PeakTimeMinusRMS();
    const double hitTimeEnd = hit->PeakTimePlusRMS();
    const double xCoord = detProp.ConvertTicksToX(hitTime, hitWireID.Plane, hitWireID.TPC, hitWireID.Cryostat);
    const double xCoordStart = detProp.ConvertTicksToX(hitTimeStart, hitWireID.Plane, hitWireID.TPC, hitWireID.Cryostat);
    const double xCoordEnd = detProp.ConvertTicksToX(hitTimeEnd, hitWireID.Plane, hitWireID.TPC, hitWireID.Cryostat);
    auto const xyz = wireReadoutGeom.Wire(hitWireID).GetCenter();
    
    pandoraHitPosition = IvysaurusUtils::ProjectIntoPandoraView(TVector3(xCoord, xyz.Y(), xyz.Z()), hitType);
    width = std::fabs(xCoordEnd - xCoordStart);
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

const TVector3 ProjectIntoPandoraView(const TVector3 &inputPosition, const PandoraView pandoraView)
{
    const float xCoord = inputPosition.X();
    const float yCoord = inputPosition.Y();
    const float zCoord = inputPosition.Z();
    return TVector3(xCoord, 0.f, pandoraView == TPC_VIEW_U ? YZToU(yCoord, zCoord) :
        pandoraView == TPC_VIEW_V ? YZToV(yCoord, zCoord) : YZToW(yCoord, zCoord));
}

//------------------------------------------------------------------------------------------------------------------------------------------

const PandoraView GetPandora2DView(const art::Ptr<recob::Hit> &hit)
{
    const geo::WireID hitWireID(hit->WireID());
    const geo::View_t hitView(hit->View());
    const geo::View_t thisPandoraView(lar_pandora::LArPandoraGeometry::GetGlobalView(hitWireID.Cryostat, hitWireID.TPC, hitView));

    if (thisPandoraView == geo::kW || thisPandoraView == geo::kY)
        return TPC_VIEW_W;
    else if (thisPandoraView == geo::kU)
        return TPC_VIEW_U;
    else if (thisPandoraView == geo::kV)
        return TPC_VIEW_V;
    else
        throw cet::exception("ivysaurus::IvysaurusUtils") << "wire view not recognised";
}

//------------------------------------------------------------------------------------------------------------------------------------------

float YZToU(const float yCoord, const float zCoord)
{
    // By definition, Pandora U/V/W wire angle is independent of iTPC/iCryo
    lar_pandora::LArPandoraDetectorType* detType(lar_pandora::detector_functions::GetDetectorType());
    const float wireAngleU(detType->WireAngleU(0, 0));
    return (zCoord * std::cos(wireAngleU)) - (yCoord * std::sin(wireAngleU));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float YZToV(const float yCoord, const float zCoord)
{
    // By definition, Pandora U/V/W wire angle is independent of iTPC/iCryo    
    lar_pandora::LArPandoraDetectorType* detType(lar_pandora::detector_functions::GetDetectorType());
    const float wireAngleV(detType->WireAngleV(0, 0));    
    return (zCoord * std::cos(wireAngleV)) - (yCoord * std::sin(wireAngleV));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float YZToW(const float yCoord, const float zCoord)
{
    // By definition, Pandora U/V/W wire angle is independent of iTPC/iCryo        
    lar_pandora::LArPandoraDetectorType* detType(lar_pandora::detector_functions::GetDetectorType());
    const float wireAngleW(detType->WireAngleW(0, 0));      
    return (zCoord * std::cos(wireAngleW)) - (yCoord * std::sin(wireAngleW));
}

//------------------------------------------------------------------------------------------------------------------------------------------

double CompletenessFromTrueParticleID(detinfo::DetectorClocksData const& clockData, const std::vector<art::Ptr<recob::Hit>> &selectedHits, 
  const std::vector<art::Ptr<recob::Hit>> &eventHits, const int trackID)
{
    int nMatchesInSelHits = 0;
    int nMatchesInAllHits = 0;

    for (art::Ptr<recob::Hit> hit : selectedHits)
    {
        const int matchedID = TruthMatchUtils::TrueParticleID(clockData, hit, 1); 
        if (matchedID == trackID) nMatchesInSelHits++;
    }

    for (art::Ptr<recob::Hit> hit : eventHits)
    {
        const int matchedID = TruthMatchUtils::TrueParticleID(clockData, hit, 1);
        if (matchedID == trackID) nMatchesInAllHits++;
    }

    const double completeness = (nMatchesInAllHits > 0) ? static_cast<double>(nMatchesInSelHits) / static_cast<double>(nMatchesInAllHits) : 0.0;

    return completeness;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double HitPurityFromTrueParticleID(detinfo::DetectorClocksData const& clockData, const std::vector<art::Ptr<recob::Hit>> &selectedHits,
    const int trackID)
{
    int nMatchesInSelHits = 0;

    for (art::Ptr<recob::Hit> hit : selectedHits)
    {
        const int matchedID = TruthMatchUtils::TrueParticleID(clockData, hit, 1);
        if (matchedID == trackID) nMatchesInSelHits++;
    }

    const double purity = (selectedHits.size() > 0) ? static_cast<double>(nMatchesInSelHits) / static_cast<double>(selectedHits.size()) : 0.0;

    return purity;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float IntegrateGaussian(const float limitA, const float limitB, const float mean, const float std, const float stepSize)
{
    const float globalIntegralMin = std::min(limitA, limitB);
    const float globalIntegralMax = std::max(limitA, limitB);
    float localIntegralMin = globalIntegralMin;
    float integralSum = 0.f;

    while (localIntegralMin < globalIntegralMax)
    {
        float localIntegralMax = localIntegralMin + stepSize;
        localIntegralMax = std::min(localIntegralMax, globalIntegralMax);
        integralSum += TrapeziumRule(localIntegralMin, localIntegralMax, mean, std);
        localIntegralMin += stepSize;
    }

    return integralSum;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TrapeziumRule(const float lowerLimit, const float upperLimit, const float mean, const float std)
{
    const float lowerExpN = (-1.f) * (lowerLimit - mean) * (lowerLimit - mean);
    const float upperExpN = (-1.f) * (upperLimit - mean) * (upperLimit - mean);
    const float expD = 2 * std * std;
    const float amp = (1.0) / (std::sqrt(2.0 * M_PI) * std);
    const float lowerEval = amp * std::exp(lowerExpN / expD);
    const float upperEval = amp * std::exp(upperExpN / expD);

    return (upperLimit - lowerLimit) * (0.5f * (lowerEval + upperEval));
}

} // namespace IvysaurusUtils

