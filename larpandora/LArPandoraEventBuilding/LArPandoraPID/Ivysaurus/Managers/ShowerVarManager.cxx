////////////////////////////////////////////////////////////////////////
/// \file    ShowerVarManager.cxx
/// \brief   A class to manage the Ivysaurus 2D shower variable input 
/// \author  Isobel Mawby - i.mawby1@lancaster.ac.uk
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <random>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larcore/Geometry/Geometry.h"

#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"

#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

#include "larpandora/LArPandoraUtils/PandoraEventUtils.h"
#include "larpandora/LArPandoraUtils/PandoraHitUtils.h"
#include "larpandora/LArPandoraUtils/PandoraPFParticleUtils.h"
#include "larpandora/LArPandoraUtils/PandoraShowerUtils.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/ShowerVarManager.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Utils/IvysaurusUtils.h"

namespace ivysaurus
{

ShowerVarManager::ShowerVars::ShowerVars() : 
    m_isNormalised(false),
    m_displacement(-1.f),
    m_DCA(-1.f),
    m_trackStubLength(-1.f),
    m_nuVertexAvSeparation(-1.f),
    m_nuVertexChargeAsymmetry(-1.f),
    m_foundConnectionPathway(0.f),
    m_initialGapSize(-1.f),
    m_largestGapSize(-1.f),
    m_pathwayLength(-1.f),
    m_pathwayScatteringAngle2D(-1.f),
    m_nShowerHits(-1.f),
    m_foundHitRatio(-1.f),
    m_scatterAngle(-1.f),
    m_openingAngle(-1.f),
    m_nuVertexEnergyAsymmetry(-1.f),
    m_nuVertexEnergyWeightedMeanRadialDistance(-1.f),
    m_showerStartEnergyAsymmetry(-1.f),
    m_showerStartMoliereRadius(-1.f),
    m_nAmbiguousViews(-1.f),
    m_unaccountedEnergy(-999.f)
{
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ShowerVarManager::ShowerVarManager(const fhicl::ParameterSet& pset) :
    m_recoModuleLabel(pset.get<std::string>("RecoModuleLabel")),
    m_showerModuleLabel(pset.get<std::string>("ShowerModuleLabel")),
    m_hitModuleLabel(pset.get<std::string>("HitModuleLabel")),
    m_displacementLimit(pset.get<float>("DisplacementLimit")),
    m_DCALimit(pset.get<float>("DCALimit")),
    m_trackStubLengthLimit(pset.get<float>("TrackStubLengthLimit")),
    m_nuVertexAvSeparationLimit(pset.get<float>("NuVertexAvSeparationLimit")),
    m_nuVertexChargeAsymmetryLimit(pset.get<float>("NuVertexChargeAsymmetryLimit")),
    m_initialGapSizeLimit(pset.get<float>("InitialGapSizeLimit")),
    m_largestGapSizeLimit(pset.get<float>("LargestGapSizeLimit")),
    m_pathwayLengthLimit(pset.get<float>("PathwayLengthLimit")),
    m_pathwayScatteringAngle2DLimit(pset.get<float>("PathwayScatteringAngle2DLimit")),
    m_nShowerHitsLimit(pset.get<float>("NShowerHitsLimit")),
    m_foundHitRatioLimit(pset.get<float>("FoundHitRatioLimit")),
    m_scatterAngleLimit(pset.get<float>("ScatterAngleLimit")),
    m_openingAngleLimit(pset.get<float>("OpeningAngleLimit")),
    m_nuVertexEnergyAsymmetryLimit(pset.get<float>("NuVertexEnergyAsymmetryLimit")),
    m_nuVertexEnergyWeightedMeanRadialDistanceLimit(pset.get<float>("NuVertexEnergyWeightedMeanRadialDistanceLimit")),
    m_showerStartEnergyAsymmetryLimit(pset.get<float>("ShowerStartEnergyAsymmetryLimit")),
    m_showerStartMoliereRadiusLimit(pset.get<float>("ShowerStartMoliereRadiusLimit")),
    m_nAmbiguousViewsLimit(pset.get<float>("NAmbiguousViewsLimit")),
    m_unaccountedEnergyLimit(pset.get<float>("UnaccountedEnergyLimit"))
{
}

/////////////////////////////////////////////////////////////

ShowerVarManager::~ShowerVarManager()
{
}

/////////////////////////////////////////////////////////////

bool ShowerVarManager::EvaluateShowerVars(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    ShowerVarManager::ShowerVars &showerVars) const
{
    if (!lar_pandora::PandoraPFParticleUtils::HasShower(pfparticle, evt, m_recoModuleLabel, m_showerModuleLabel))
        return false;

    const art::Ptr<recob::Shower> shower = lar_pandora::PandoraPFParticleUtils::GetShower(pfparticle, evt, m_recoModuleLabel, m_showerModuleLabel);

    FillDisplacement(evt, shower, showerVars);
    FillConnectionPathwayVars(evt, pfparticle, showerVars);
    FillTrackStub(evt, shower, showerVars);
    FillNuVertexAvSeparation(evt, pfparticle, showerVars);
    FillNuVertexChargeAsymmetry(evt, shower, showerVars);

    return true;
}

/////////////////////////////////////////////////////////////

void ShowerVarManager::FillDisplacement(const art::Event &evt, const art::Ptr<recob::Shower> &shower, 
    ShowerVarManager::ShowerVars &showerVars) const
{
    if (!lar_pandora::PandoraEventUtils::GetNeutrino(evt, m_recoModuleLabel))
        return;

    const art::Ptr<recob::PFParticle> &nuPFP = lar_pandora::PandoraEventUtils::GetNeutrino(evt, m_recoModuleLabel);
    const art::Ptr<recob::Vertex> &nuVertex = lar_pandora::PandoraPFParticleUtils::GetVertex(nuPFP, evt, m_recoModuleLabel);
    const TVector3 nuVertexPosition = TVector3(nuVertex->position().X(), nuVertex->position().Y(), nuVertex->position().Z());
    const TVector3 showerStart = TVector3(shower->ShowerStart().X(), shower->ShowerStart().Y(), shower->ShowerStart().Z());
    const float displacement = (nuVertexPosition - showerStart).Mag();

    // DCA
    const double alpha = std::fabs((shower->ShowerStart() - nuVertexPosition).Dot(shower->Direction()));
    const TVector3 r = shower->ShowerStart() - (alpha * shower->Direction());
    const float dca = (r - nuVertexPosition).Mag();

    showerVars.SetDisplacement(displacement);
    showerVars.SetDCA(dca);
}

/////////////////////////////////////////////////////////////

void ShowerVarManager::FillConnectionPathwayVars(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    ShowerVarManager::ShowerVars &showerVars) const
{
    const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata = lar_pandora::PandoraPFParticleUtils::GetMetadata(pfparticle, evt, m_recoModuleLabel);
    const larpandoraobj::PFParticleMetadata::PropertiesMap &propertyMap = metadata->GetPropertiesMap();

    if (propertyMap.find("FoundConnectionPathway") != propertyMap.end()) showerVars.SetFoundConnectionPathway(propertyMap.at("FoundConnectionPathway") ? 1.f : 0.f);
    if (propertyMap.find("MaxInitialGapSize") != propertyMap.end()) showerVars.SetInitialGapSize(propertyMap.at("MaxInitialGapSize"));
    if (propertyMap.find("MinLargestProjectedGapSize") != propertyMap.end()) showerVars.SetLargestGapSize(propertyMap.at("MinLargestProjectedGapSize"));
    if (propertyMap.find("PathwayLengthMin") != propertyMap.end()) showerVars.SetPathwayLength(propertyMap.at("PathwayLengthMin"));
    if (propertyMap.find("MaxShowerStartPathwayScatteringAngle2D") != propertyMap.end()) showerVars.SetPathwayScatteringAngle2D(propertyMap.at("MaxShowerStartPathwayScatteringAngle2D"));
    if (propertyMap.find("MaxNPostShowerStartHits") != propertyMap.end()) showerVars.SetNShowerHits(propertyMap.at("MaxNPostShowerStartHits"));
    if (propertyMap.find("MaxFoundHitRatio") != propertyMap.end()) showerVars.SetFoundHitRatio(propertyMap.at("MaxFoundHitRatio"));
    if (propertyMap.find("MaxPostShowerStartScatterAngle") != propertyMap.end()) showerVars.SetScatterAngle(propertyMap.at("MaxPostShowerStartScatterAngle"));
    if (propertyMap.find("MaxPostShowerStartOpeningAngle") != propertyMap.end()) showerVars.SetOpeningAngle(propertyMap.at("MaxPostShowerStartOpeningAngle"));
    if (propertyMap.find("MaxPostShowerStartNuVertexEnergyAsymmetry") != propertyMap.end()) showerVars.SetNuVertexEnergyAsymmetry(propertyMap.at("MaxPostShowerStartNuVertexEnergyAsymmetry"));
    if (propertyMap.find("MaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance") != propertyMap.end()) 
        showerVars.SetNuVertexEnergyWeightedMeanRadialDistance(propertyMap.at("MaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance"));
    if (propertyMap.find("MaxPostShowerStartShowerStartEnergyAsymmetry") != propertyMap.end()) showerVars.SetShowerStartEnergyAsymmetry(propertyMap.at("MaxPostShowerStartShowerStartEnergyAsymmetry"));
    if (propertyMap.find("MinPostShowerStartShowerStartMoliereRadius") != propertyMap.end()) showerVars.SetShowerStartMoliereRadius(propertyMap.at("MinPostShowerStartShowerStartMoliereRadius"));
    if (propertyMap.find("NViewsWithAmbiguousHits") != propertyMap.end()) showerVars.SetNAmbiguousViews(propertyMap.at("NViewsWithAmbiguousHits"));
    if (propertyMap.find("AmbiguousHitMaxUnaccountedEnergy") != propertyMap.end()) showerVars.SetUnaccountedEnergy(propertyMap.at("AmbiguousHitMaxUnaccountedEnergy"));
}

/////////////////////////////////////////////////////////////

void ShowerVarManager::FillTrackStub(const art::Event &evt, const art::Ptr<recob::Shower> shower, 
    ShowerVarManager::ShowerVars &showerVars) const
{
    art::Handle< std::vector<recob::Shower> > showerListHandle;
    evt.getByLabel(m_showerModuleLabel, showerListHandle);

    art::FindManyP<recob::Track> initialTrackAssn(showerListHandle, evt, m_showerModuleLabel);
    std::vector<art::Ptr<recob::Track>> initialTrackStub = initialTrackAssn.at(shower.key());

    if (initialTrackStub.empty())
        return;

    const art::Ptr<recob::Track> &trackStub = initialTrackStub.front();

    const float trackStubLength = std::sqrt((trackStub->Start() - trackStub->End()).Mag2());
    showerVars.SetTrackStubLength(trackStubLength);
}

/////////////////////////////////////////////////////////////

void ShowerVarManager::FillNuVertexAvSeparation(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    ShowerVarManager::ShowerVars &showerVars) const
{
    if (!lar_pandora::PandoraPFParticleUtils::IsShower(pfparticle, evt, m_recoModuleLabel, m_showerModuleLabel))
        return;

    const art::Ptr<recob::Shower> shower = lar_pandora::PandoraPFParticleUtils::GetShower(pfparticle, evt, m_recoModuleLabel, m_showerModuleLabel);

    if (!lar_pandora::PandoraEventUtils::GetNeutrino(evt, m_recoModuleLabel))
        return;

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

    const art::Ptr<recob::PFParticle> &nuPFP = lar_pandora::PandoraEventUtils::GetNeutrino(evt, m_recoModuleLabel);
    const art::Ptr<recob::Vertex> &nuVertex = lar_pandora::PandoraPFParticleUtils::GetVertex(nuPFP, evt, m_recoModuleLabel);
    const TVector3 nuVertexPosition = TVector3(nuVertex->position().X(), nuVertex->position().Y(), nuVertex->position().Z());
    const TVector3 showerStart = TVector3(shower->ShowerStart().X(), shower->ShowerStart().Y(), shower->ShowerStart().Z());
    const TVector3 displacement = (showerStart - nuVertexPosition).Unit();

    // Calc charge weighted axis-hit separation
    float totalCharge = 0.f;
    float numeratorSum = 0.f;

    const std::vector<art::Ptr<recob::Hit>> collectionHits = lar_pandora::PandoraPFParticleUtils::GetViewHits(pfparticle, evt, m_recoModuleLabel, 2);

    for (const art::Ptr<recob::Hit> &hit : collectionHits)
    {
        const std::vector<art::Ptr<recob::SpacePoint>> &spacepoints = lar_pandora::PandoraHitUtils::GetSpacePoints(hit, evt, m_hitModuleLabel, m_recoModuleLabel);

        if (spacepoints.empty())
            continue;

        const art::Ptr<recob::SpacePoint> spacepoint = spacepoints.front();

        const TVector3 spacepointPos = TVector3(spacepoint->position().X(), spacepoint->position().Y(), spacepoint->position().Z()) - nuVertexPosition;
        const float transverse = std::sqrt(displacement.Cross(spacepointPos).Mag2());
        const double charge = lar_pandora::PandoraHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, {hit});

        totalCharge += charge;
        numeratorSum += (transverse * charge);
    }

    const float nuVertexAvSeparation = (totalCharge < std::numeric_limits<float>::epsilon()) ? 0.f : (numeratorSum / totalCharge);
    showerVars.SetNuVertexAvSeparation(nuVertexAvSeparation);
}

/////////////////////////////////////////////////////////////

void ShowerVarManager::FillNuVertexChargeAsymmetry(const art::Event &evt, const art::Ptr<recob::Shower> &shower, 
    ShowerVarManager::ShowerVars &showerVars) const
{
    if (!lar_pandora::PandoraEventUtils::GetNeutrino(evt, m_recoModuleLabel))
        return;

    const art::Ptr<recob::PFParticle> &nuPFP = lar_pandora::PandoraEventUtils::GetNeutrino(evt, m_recoModuleLabel);
    const art::Ptr<recob::Vertex> &nuVertex = lar_pandora::PandoraPFParticleUtils::GetVertex(nuPFP, evt, m_recoModuleLabel);
    const TVector3 nuVertexPosition = TVector3(nuVertex->position().X(), nuVertex->position().Y(), nuVertex->position().Z());

    // Get initial direction of the shower from the track stub
    const TVector3 showerDirection = TVector3(shower->Direction().X(), shower->Direction().Y(), shower->Direction().Z());

    const std::vector<art::Ptr<recob::Hit>> &hits = lar_pandora::PandoraShowerUtils::GetHits(shower, evt, m_showerModuleLabel);

    std::vector<art::Ptr<recob::Hit>> hitsU, hitsV, hitsW;

    for (const art::Ptr<recob::Hit> &hit : hits)
    {
        // Make sure 2D hit has an associated space point
        std::vector<art::Ptr<recob::SpacePoint>> spacePoints = lar_pandora::PandoraHitUtils::GetSpacePoints(hit, evt, m_hitModuleLabel, m_recoModuleLabel);

        if (spacePoints.empty())
            continue;

        const IvysaurusUtils::PandoraView pandoraView = IvysaurusUtils::GetPandora2DView(hit);

        std::vector<art::Ptr<recob::Hit>> &hitVector = pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_U ? hitsU :
            pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_V ? hitsV : hitsW;

        hitVector.push_back(hit);
    }

    const float nuVertexChargeAsymmetryU = GetViewNuVertexChargeAsymmetry(evt, nuVertexPosition, showerDirection, hitsU, IvysaurusUtils::PandoraView::TPC_VIEW_U);
    const float nuVertexChargeAsymmetryV = GetViewNuVertexChargeAsymmetry(evt, nuVertexPosition, showerDirection, hitsV, IvysaurusUtils::PandoraView::TPC_VIEW_V);
    const float nuVertexChargeAsymmetryW = GetViewNuVertexChargeAsymmetry(evt, nuVertexPosition, showerDirection, hitsW, IvysaurusUtils::PandoraView::TPC_VIEW_W);
    const float maxNuVertexChargeAsymmetry = std::max(std::max(nuVertexChargeAsymmetryU, nuVertexChargeAsymmetryV), nuVertexChargeAsymmetryW);

    showerVars.SetNuVertexChargeAsymmetry(maxNuVertexChargeAsymmetry);
}

/////////////////////////////////////////////////////////////

float ShowerVarManager::GetViewNuVertexChargeAsymmetry(const art::Event &evt, const TVector3 &nuVertexPosition, const TVector3 &showerDirection, 
    const std::vector<art::Ptr<recob::Hit>> &viewHits, const IvysaurusUtils::PandoraView &pandoraView) const
{
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

    const TVector3 viewNuVertex = IvysaurusUtils::ProjectIntoPandoraView(nuVertexPosition, pandoraView);
    const TVector3 point2 = nuVertexPosition + (showerDirection * 10.f);
    const TVector3 viewPoint2 = IvysaurusUtils::ProjectIntoPandoraView(point2, pandoraView);
    const TVector3 centralAxis = (viewPoint2 - viewNuVertex).Unit();
    const TVector3 yAxis = TVector3(0.f, 1.f, 0.f);
    const TVector3 orthAxis = centralAxis.Cross(yAxis).Unit();

    float chargeAsymmetry = 0.f;
    float totalCharge = 0.f;

    for (const art::Ptr<recob::Hit> &viewHit : viewHits)
    {
        float width = 0.f;
        TVector3 hitPosition = TVector3(0.f, 0.f, 0.f);
        IvysaurusUtils::ObtainPandoraHitPositionAndWidth(evt, viewHit, pandoraView, hitPosition, width);
        const float l = orthAxis.Dot(hitPosition - viewNuVertex);
        float charge = std::fabs(lar_pandora::PandoraHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, {viewHit}));

        totalCharge += charge;

        charge *= (l < 0.f) ? -1.0 : 1.0;
        chargeAsymmetry += charge;
    }

    chargeAsymmetry = totalCharge < std::numeric_limits<float>::epsilon() ? 0.f : (chargeAsymmetry / totalCharge);

    return std::fabs(chargeAsymmetry);
}

/////////////////////////////////////////////////////////////

void ShowerVarManager::Reset(ShowerVarManager::ShowerVars &showerVars) const
{
    showerVars.SetIsNormalised(false);
    showerVars.SetDisplacement(-1.f);
    showerVars.SetDCA(-1.f);
    showerVars.SetTrackStubLength(-1.f);
    showerVars.SetNuVertexAvSeparation(-1.f);
    showerVars.SetNuVertexChargeAsymmetry(-1.f);
    showerVars.SetInitialGapSize(-1.f);
    showerVars.SetLargestGapSize(-1.f);
    showerVars.SetPathwayLength(-1.f);
    showerVars.SetPathwayScatteringAngle2D(-1.f);
    showerVars.SetNShowerHits(-1.f);
    showerVars.SetFoundHitRatio(-1.f);
    showerVars.SetScatterAngle(-1.f);
    showerVars.SetOpeningAngle(-1.f);
    showerVars.SetNuVertexEnergyAsymmetry(-1.f);
    showerVars.SetNuVertexEnergyWeightedMeanRadialDistance(-1.f);
    showerVars.SetShowerStartEnergyAsymmetry(-1.f);
    showerVars.SetShowerStartMoliereRadius(-1.f);
    showerVars.SetNAmbiguousViews(-1.f);
    showerVars.SetUnaccountedEnergy(-999.f);
}

/////////////////////////////////////////////////////////////

void ShowerVarManager::NormaliseShowerVars(ShowerVarManager::ShowerVars &showerVars) const
{
    if (showerVars.GetIsNormalised())
        return;

    // Displacement
    if (showerVars.GetDisplacement() > m_displacementLimit)
        showerVars.SetDisplacement(m_displacementLimit);

    if (showerVars.GetDisplacement() < ((-1.0) * std::numeric_limits<float>::epsilon()))
        showerVars.SetDisplacement((-1.0) * m_displacementLimit);

    showerVars.SetDisplacement(showerVars.GetDisplacement() / m_displacementLimit);

    // DCA
    if (showerVars.GetDCA() > m_DCALimit)
        showerVars.SetDCA(m_DCALimit);

    if (showerVars.GetDCA() < ((-1.0) * std::numeric_limits<float>::epsilon()))
        showerVars.SetDCA((-1.0) * m_DCALimit);

    showerVars.SetDCA(showerVars.GetDCA() / m_DCALimit);

    // TrackStubLength
    if (showerVars.GetTrackStubLength() > m_trackStubLengthLimit)
        showerVars.SetTrackStubLength(m_trackStubLengthLimit);

    if (showerVars.GetTrackStubLength() < ((-1.0) * std::numeric_limits<float>::epsilon()))
        showerVars.SetTrackStubLength((-1.0) * m_trackStubLengthLimit);

    showerVars.SetTrackStubLength(showerVars.GetTrackStubLength() / m_trackStubLengthLimit);

    // NuVertexAvSep
    if (showerVars.GetNuVertexAvSeparation() > m_nuVertexAvSeparationLimit)
        showerVars.SetNuVertexAvSeparation(m_nuVertexAvSeparationLimit);

    if (showerVars.GetNuVertexAvSeparation() < ((-1.0) * std::numeric_limits<float>::epsilon()))
        showerVars.SetNuVertexAvSeparation((-1.0) * m_nuVertexAvSeparationLimit);

    showerVars.SetNuVertexAvSeparation(showerVars.GetNuVertexAvSeparation() / m_nuVertexAvSeparationLimit);

    // NuVertexChargeAsymmetry
    if (showerVars.GetNuVertexChargeAsymmetry() > m_nuVertexChargeAsymmetryLimit)
        showerVars.SetNuVertexChargeAsymmetry(m_nuVertexChargeAsymmetryLimit);

    if (showerVars.GetNuVertexChargeAsymmetry() < ((-1.0) * std::numeric_limits<float>::epsilon()))
        showerVars.SetNuVertexChargeAsymmetry((-1.0) * m_nuVertexChargeAsymmetryLimit);

    showerVars.SetNuVertexChargeAsymmetry(showerVars.GetNuVertexChargeAsymmetry() / m_nuVertexChargeAsymmetryLimit);

    // InitialGapSize
    if (showerVars.GetInitialGapSize() > m_initialGapSizeLimit)
        showerVars.SetInitialGapSize(m_initialGapSizeLimit);

    if (showerVars.GetInitialGapSize() < ((-1.0) * std::numeric_limits<float>::epsilon()))
        showerVars.SetInitialGapSize((-1.0) * m_initialGapSizeLimit);

    showerVars.SetInitialGapSize(showerVars.GetInitialGapSize() / m_initialGapSizeLimit);

    // LargestGapSize
    if (showerVars.GetLargestGapSize() > m_largestGapSizeLimit)
        showerVars.SetLargestGapSize(m_largestGapSizeLimit);

    if (showerVars.GetLargestGapSize() < ((-1.0) * std::numeric_limits<float>::epsilon()))
        showerVars.SetLargestGapSize((-1.0) * m_largestGapSizeLimit);

    showerVars.SetLargestGapSize(showerVars.GetLargestGapSize() / m_largestGapSizeLimit);

    // PathwayLength
    if (showerVars.GetPathwayLength() > m_pathwayLengthLimit)
        showerVars.SetPathwayLength(m_pathwayLengthLimit);

    if (showerVars.GetPathwayLength() < ((-1.0) * std::numeric_limits<float>::epsilon()))
        showerVars.SetPathwayLength((-1.0) * m_pathwayLengthLimit);

    showerVars.SetPathwayLength(showerVars.GetPathwayLength() / m_pathwayLengthLimit);

    // PathwayScatteringAngle2D
    if (showerVars.GetPathwayScatteringAngle2D() > m_pathwayScatteringAngle2DLimit)
        showerVars.SetPathwayScatteringAngle2D(m_pathwayScatteringAngle2DLimit);

    if (showerVars.GetPathwayScatteringAngle2D() < ((-1.0) * std::numeric_limits<float>::epsilon()))
        showerVars.SetPathwayScatteringAngle2D((-1.0) * m_pathwayScatteringAngle2DLimit);

    showerVars.SetPathwayScatteringAngle2D(showerVars.GetPathwayScatteringAngle2D() / m_pathwayScatteringAngle2DLimit);

    // NShowerHits
    if (showerVars.GetNShowerHits() > m_nShowerHitsLimit)
        showerVars.SetNShowerHits(m_nShowerHitsLimit);

    if (showerVars.GetNShowerHits() < ((-1.0) * std::numeric_limits<float>::epsilon()))
        showerVars.SetNShowerHits((-1.0) * m_nShowerHitsLimit);

    showerVars.SetNShowerHits(showerVars.GetNShowerHits() / m_nShowerHitsLimit);

    // FoundHitRatio
    if (showerVars.GetFoundHitRatio() > m_foundHitRatioLimit)
        showerVars.SetFoundHitRatio(m_foundHitRatioLimit);

    if (showerVars.GetFoundHitRatio() < ((-1.0) * std::numeric_limits<float>::epsilon()))
        showerVars.SetFoundHitRatio((-1.0) * m_foundHitRatioLimit);

    showerVars.SetFoundHitRatio(showerVars.GetFoundHitRatio() / m_foundHitRatioLimit);

    // ScatterAngle
    if (showerVars.GetScatterAngle() > m_scatterAngleLimit)
        showerVars.SetScatterAngle(m_scatterAngleLimit);

    if (showerVars.GetScatterAngle() < ((-1.0) * std::numeric_limits<float>::epsilon()))
        showerVars.SetScatterAngle((-1.0) * m_scatterAngleLimit);

    showerVars.SetScatterAngle(showerVars.GetScatterAngle() / m_scatterAngleLimit);

    // OpeningAngle
    if (showerVars.GetOpeningAngle() > m_openingAngleLimit)
        showerVars.SetOpeningAngle(m_openingAngleLimit);

    if (showerVars.GetOpeningAngle() < ((-1.0) * std::numeric_limits<float>::epsilon()))
        showerVars.SetOpeningAngle((-1.0) * m_openingAngleLimit);

    showerVars.SetOpeningAngle(showerVars.GetOpeningAngle() / m_openingAngleLimit);

    // NuVertexEnergyAsymmetry
    if (showerVars.GetNuVertexEnergyAsymmetry() > m_nuVertexEnergyAsymmetryLimit)
        showerVars.SetNuVertexEnergyAsymmetry(m_nuVertexEnergyAsymmetryLimit);

    if (showerVars.GetNuVertexEnergyAsymmetry() < ((-1.0) * std::numeric_limits<float>::epsilon()))
        showerVars.SetNuVertexEnergyAsymmetry((-1.0) * m_nuVertexEnergyAsymmetryLimit);

    showerVars.SetNuVertexEnergyAsymmetry(showerVars.GetNuVertexEnergyAsymmetry() / m_nuVertexEnergyAsymmetryLimit);

    // NuVertexEnergyWeightedMeanRadialDistance
    if (showerVars.GetNuVertexEnergyWeightedMeanRadialDistance() > m_nuVertexEnergyWeightedMeanRadialDistanceLimit)
        showerVars.SetNuVertexEnergyWeightedMeanRadialDistance(m_nuVertexEnergyWeightedMeanRadialDistanceLimit);

    if (showerVars.GetNuVertexEnergyWeightedMeanRadialDistance() < ((-1.0) * std::numeric_limits<float>::epsilon()))
        showerVars.SetNuVertexEnergyWeightedMeanRadialDistance((-1.0) * m_nuVertexEnergyWeightedMeanRadialDistanceLimit);

    showerVars.SetNuVertexEnergyWeightedMeanRadialDistance(showerVars.GetNuVertexEnergyWeightedMeanRadialDistance() / m_nuVertexEnergyWeightedMeanRadialDistanceLimit);

    // ShowerStartEnergyAsymmetry
    if (showerVars.GetShowerStartEnergyAsymmetry() > m_showerStartEnergyAsymmetryLimit)
        showerVars.SetShowerStartEnergyAsymmetry(m_showerStartEnergyAsymmetryLimit);

    if (showerVars.GetShowerStartEnergyAsymmetry() < ((-1.0) * std::numeric_limits<float>::epsilon()))
        showerVars.SetShowerStartEnergyAsymmetry((-1.0) * m_showerStartEnergyAsymmetryLimit);

    showerVars.SetShowerStartEnergyAsymmetry(showerVars.GetShowerStartEnergyAsymmetry() / m_showerStartEnergyAsymmetryLimit);

    // ShowerStartMoliereRadius
    if (showerVars.GetShowerStartMoliereRadius() > m_showerStartMoliereRadiusLimit)
        showerVars.SetShowerStartMoliereRadius(m_showerStartMoliereRadiusLimit);

    if (showerVars.GetShowerStartMoliereRadius() < ((-1.0) * std::numeric_limits<float>::epsilon()))
        showerVars.SetShowerStartMoliereRadius((-1.0) * m_showerStartMoliereRadiusLimit);

    showerVars.SetShowerStartMoliereRadius(showerVars.GetShowerStartMoliereRadius() / m_showerStartMoliereRadiusLimit);

    // NAmbiguousViews
    if (showerVars.GetNAmbiguousViews() > m_nAmbiguousViewsLimit)
        showerVars.SetNAmbiguousViews(m_nAmbiguousViewsLimit);

    if (showerVars.GetNAmbiguousViews() < ((-1.0) * std::numeric_limits<float>::epsilon()))
        showerVars.SetNAmbiguousViews((-1.0) * m_nAmbiguousViewsLimit);

    showerVars.SetNAmbiguousViews(showerVars.GetNAmbiguousViews() / m_nAmbiguousViewsLimit);

    // UnaccountedEnergy
    if (showerVars.GetUnaccountedEnergy() > m_unaccountedEnergyLimit)
        showerVars.SetUnaccountedEnergy(m_unaccountedEnergyLimit);

    if (showerVars.GetUnaccountedEnergy() < ((-1.0) * m_unaccountedEnergyLimit))
        showerVars.SetUnaccountedEnergy((-1.0) * m_unaccountedEnergyLimit);

    showerVars.SetUnaccountedEnergy(showerVars.GetUnaccountedEnergy() / m_unaccountedEnergyLimit);

    // Set normalised
    showerVars.SetIsNormalised(true);
}

} //namespace ivysaurus
