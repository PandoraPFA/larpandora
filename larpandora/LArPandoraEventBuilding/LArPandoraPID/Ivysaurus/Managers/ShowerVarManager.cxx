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
    m_displacement(std::make_pair(-1.f, false)),
    m_DCA(std::make_pair(-1.f, false)),
    m_trackStubLength(std::make_pair(-1.f, false)),
    m_nuVertexAvSeparation(std::make_pair(-1.f, false)),
    m_nuVertexChargeAsymmetry(std::make_pair(-1.f, false))
{
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ShowerVarManager::ShowerVarManager(const fhicl::ParameterSet& pset) :
    m_recoModuleLabel(pset.get<std::string>("RecoModuleLabel")),
    m_trackModuleLabel(pset.get<std::string>("TrackModuleLabel")),    
    m_showerModuleLabel(pset.get<std::string>("ShowerModuleLabel")),
    m_hitModuleLabel(pset.get<std::string>("HitModuleLabel")),
    m_displacementMean(pset.get<float>("DisplacementMean")),
    m_displacementStd(pset.get<float>("DisplacementStd")),
    m_DCAMean(pset.get<float>("DCAMean")),
    m_DCAStd(pset.get<float>("DCAStd")),
    m_trackStubLengthMean(pset.get<float>("TrackStubLengthMean")),
    m_trackStubLengthStd(pset.get<float>("TrackStubLengthStd")),
    m_nuVertexAvSeparationMean(pset.get<float>("NuVertexAvSeparationMean")),
    m_nuVertexAvSeparationStd(pset.get<float>("NuVertexAvSeparationStd")), 
    m_nuVertexChargeAsymmetryMean(pset.get<float>("NuVertexChargeAsymmetryMean")),
    m_nuVertexChargeAsymmetryStd(pset.get<float>("NuVertexChargeAsymmetryStd"))    
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

    // Network not primary particle specific
    TVector3 parentEndpoint(0.f, 0.f, 0.f);
    if (!this->GetParentEndpoint(evt, pfparticle, parentEndpoint))
        return false;

    const art::Ptr<recob::Shower> shower = lar_pandora::PandoraPFParticleUtils::GetShower(pfparticle, evt, m_recoModuleLabel, m_showerModuleLabel);
    FillDisplacement(evt, parentEndpoint, shower, showerVars);
    FillTrackStub(evt, shower, showerVars);
    FillAvSeparation(evt, parentEndpoint, pfparticle, shower, showerVars);
    FillChargeAsymmetry(evt, parentEndpoint, shower, showerVars);

    return true;
}

/////////////////////////////////////////////////////////////

bool ShowerVarManager::GetParentEndpoint(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, TVector3 &parentEnd) const
{
    // Find the parent
    bool found(false);
    art::Ptr<recob::PFParticle> parentPFP;
    unsigned int parentSelf(pfparticle->Parent());
    const std::vector<art::Ptr<recob::PFParticle>> &pfps = lar_pandora::PandoraEventUtils::GetPFParticles(evt, m_recoModuleLabel);
    
    for (art::Ptr<recob::PFParticle> pfp : pfps)
    {
        if (pfp->Self() == parentSelf)
        {
            parentPFP = pfp;
            found = true;
            break;
        }
    }

    if (!found)
        return false;

    // Get endpoint
    if ((parentPFP->PdgCode() == 12) || (parentPFP->PdgCode() == 14) || (parentPFP->PdgCode() == 16))
    {
        const art::Ptr<recob::Vertex> &parentVertex = lar_pandora::PandoraPFParticleUtils::GetVertex(parentPFP, evt, m_recoModuleLabel);
        const TVector3 parentVertexPosition = TVector3(parentVertex->position().X(), parentVertex->position().Y(), parentVertex->position().Z());
        parentEnd = parentVertexPosition;
    }
    else
    {
        if (!lar_pandora::PandoraPFParticleUtils::HasTrack(parentPFP, evt, m_recoModuleLabel, m_trackModuleLabel))
            return false;

        const art::Ptr<recob::Track> track = lar_pandora::PandoraPFParticleUtils::GetTrack(parentPFP, evt, m_recoModuleLabel, m_trackModuleLabel);
        parentEnd = TVector3(track->End().X(), track->End().Y(), track->End().Z());
    }

    return true;
}

/////////////////////////////////////////////////////////////

void ShowerVarManager::FillDisplacement(const art::Event &evt, const TVector3 &parentEndpoint, const art::Ptr<recob::Shower> &shower, 
    ShowerVarManager::ShowerVars &showerVars) const
{
    const TVector3 showerStart = TVector3(shower->ShowerStart().X(), shower->ShowerStart().Y(), shower->ShowerStart().Z());
    const float displacement = (parentEndpoint - showerStart).Mag();

    // DCA
    const double alpha = std::fabs((shower->ShowerStart() - parentEndpoint).Dot(shower->Direction()));
    const TVector3 r = shower->ShowerStart() - (alpha * shower->Direction());
    const float dca = (r - parentEndpoint).Mag();

    showerVars.SetDisplacement(displacement);
    showerVars.SetDCA(dca);
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

void ShowerVarManager::FillAvSeparation(const art::Event &evt, const TVector3 &parentEndpoint, const art::Ptr<recob::PFParticle> &pfparticle, 
    const art::Ptr<recob::Shower> &shower, ShowerVarManager::ShowerVars &showerVars) const
{
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

    const TVector3 showerStart = TVector3(shower->ShowerStart().X(), shower->ShowerStart().Y(), shower->ShowerStart().Z());
    const TVector3 displacement = (showerStart - parentEndpoint).Unit();

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

        const TVector3 spacepointPos = TVector3(spacepoint->position().X(), spacepoint->position().Y(), spacepoint->position().Z()) - parentEndpoint;
        const float transverse = std::sqrt(displacement.Cross(spacepointPos).Mag2());
        const double charge = lar_pandora::PandoraHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, {hit});

        totalCharge += charge;
        numeratorSum += (transverse * charge);
    }

    const float nuVertexAvSeparation = (totalCharge < std::numeric_limits<float>::epsilon()) ? 0.f : (numeratorSum / totalCharge);
    showerVars.SetNuVertexAvSeparation(nuVertexAvSeparation);
}

/////////////////////////////////////////////////////////////

void ShowerVarManager::FillChargeAsymmetry(const art::Event &evt, const TVector3 &parentEndpoint, const art::Ptr<recob::Shower> &shower, 
    ShowerVarManager::ShowerVars &showerVars) const
{
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

    const float chargeAsymmetryU = GetViewChargeAsymmetry(evt, parentEndpoint, showerDirection, hitsU, IvysaurusUtils::PandoraView::TPC_VIEW_U);
    const float chargeAsymmetryV = GetViewChargeAsymmetry(evt, parentEndpoint, showerDirection, hitsV, IvysaurusUtils::PandoraView::TPC_VIEW_V);
    const float chargeAsymmetryW = GetViewChargeAsymmetry(evt, parentEndpoint, showerDirection, hitsW, IvysaurusUtils::PandoraView::TPC_VIEW_W);
    const float maxChargeAsymmetry = std::max(std::max(chargeAsymmetryU, chargeAsymmetryV), chargeAsymmetryW);

    showerVars.SetNuVertexChargeAsymmetry(maxChargeAsymmetry);
}

/////////////////////////////////////////////////////////////

float ShowerVarManager::GetViewChargeAsymmetry(const art::Event &evt, const TVector3 &parentEndpoint, const TVector3 &showerDirection, 
    const std::vector<art::Ptr<recob::Hit>> &viewHits, const IvysaurusUtils::PandoraView &pandoraView) const
{
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

    const TVector3 viewEndpoint = IvysaurusUtils::ProjectIntoPandoraView(parentEndpoint, pandoraView);
    const TVector3 point2 = parentEndpoint + (showerDirection * 10.f);
    const TVector3 viewPoint2 = IvysaurusUtils::ProjectIntoPandoraView(point2, pandoraView);
    const TVector3 centralAxis = (viewPoint2 - viewEndpoint).Unit();
    const TVector3 yAxis = TVector3(0.f, 1.f, 0.f);
    const TVector3 orthAxis = centralAxis.Cross(yAxis).Unit();

    float chargeAsymmetry = 0.f;
    float totalCharge = 0.f;

    for (const art::Ptr<recob::Hit> &viewHit : viewHits)
    {
        float width = 0.f;
        TVector3 hitPosition = TVector3(0.f, 0.f, 0.f);
        IvysaurusUtils::ObtainPandoraHitPositionAndWidth(evt, viewHit, pandoraView, hitPosition, width);
        const float l = orthAxis.Dot(hitPosition - viewEndpoint);
        float charge = std::fabs(lar_pandora::PandoraHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, {viewHit}));

        totalCharge += charge;

        charge *= (l < 0.f) ? -1.0 : 1.0;
        chargeAsymmetry += charge;
    }

    chargeAsymmetry = totalCharge < std::numeric_limits<float>::epsilon() ? 0.f : (chargeAsymmetry / totalCharge);

    return std::fabs(chargeAsymmetry);
}

/////////////////////////////////////////////////////////////

void ShowerVarManager::NormaliseShowerVars(ShowerVarManager::ShowerVars &showerVars) const
{
    if (showerVars.GetIsNormalised())
        throw cet::exception("ivysaur::TrackVarManager::NormaliseShowerVars: ") << "shower vars are already normalised!";  

    if (showerVars.GetDisplacement().second)
        showerVars.SetDisplacement(this->NormaliseShowerVar(showerVars.GetDisplacement(), m_displacementMean, m_displacementStd));

    if (showerVars.GetDCA().second)
        showerVars.SetDCA(this->NormaliseShowerVar(showerVars.GetDCA(), m_DCAMean, m_DCAStd));

    if (showerVars.GetTrackStubLength().second)
        showerVars.SetTrackStubLength(this->NormaliseShowerVar(showerVars.GetTrackStubLength(), m_trackStubLengthMean, m_trackStubLengthStd));

    if (showerVars.GetNuVertexAvSeparation().second)
        showerVars.SetNuVertexAvSeparation(this->NormaliseShowerVar(showerVars.GetNuVertexAvSeparation(), m_nuVertexAvSeparationMean, m_nuVertexAvSeparationStd));

    if (showerVars.GetNuVertexChargeAsymmetry().second)
        showerVars.SetNuVertexChargeAsymmetry(this->NormaliseShowerVar(showerVars.GetNuVertexChargeAsymmetry(), m_nuVertexChargeAsymmetryMean, m_nuVertexChargeAsymmetryStd));

    showerVars.SetIsNormalised(true);
}

/////////////////////////////////////////////////////////////    
    
float ShowerVarManager::NormaliseShowerVar(const std::pair<float, bool> &inputShowerVar, const float mean, const float std) const
{
    // Is still on default value...
    if (!inputShowerVar.second)
        throw cet::exception("ivysaur::ShowerVarManager::NormaliseShowerVar: ") << "can't normalise a variables that wasn't set";
    
    const float normalised((inputShowerVar.first - mean) / std);
    return normalised;
}

/////////////////////////////////////////////////////////////    

} //namespace ivysaurus
