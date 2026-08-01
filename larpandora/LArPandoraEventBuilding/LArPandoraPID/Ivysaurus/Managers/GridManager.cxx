/**
 *  @file  larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/GridManager.h
 *
 *  @brief A class to manage the Ivysaurus 2D grid input 
 *
 */
// ART
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
// LArSoft
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/GridManager.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Utils/IvysaurusUtils.h"
#include "larpandora/LArPandoraUtils/PandoraEventUtils.h"
#include "larpandora/LArPandoraUtils/PandoraPFParticleUtils.h"
#include "larpandora/LArPandoraUtils/PandoraSpacePointUtils.h"
#include "larpandora/LArPandoraUtils/PandoraHitUtils.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
// Pandora
#include "Pandora/PdgTable.h"

namespace ivysaurus
{

GridManager::Grid::Grid(const TVector3 origin, const float driftSpan, const float wireSpan, 
    const unsigned int dimensions,const unsigned int nSigmaConsidered, const float integralStep) : 
        m_axisDimensions(dimensions),
        m_nSigmaConsidered(nSigmaConsidered),
        m_integralStep(integralStep)
{
    if (m_axisDimensions == 0)
        throw cet::exception("ivysaurus::Grid") << "axisDimensions can't be zero!";

    m_gridValues = GridManager::GridValues(m_axisDimensions, std::vector<std::pair<float, bool>>(m_axisDimensions, std::make_pair(0.f, false)));    

    const float driftInterval = driftSpan / m_axisDimensions;

    for (unsigned int i = 0; i <= m_axisDimensions; ++i)
        m_driftBoundaries.push_back(origin.X() + (i * driftInterval));

    const float wireInterval = wireSpan / m_axisDimensions;

    for (unsigned int i = 0; i <= m_axisDimensions; ++i)
        m_wireBoundaries.push_back(origin.Z() + (i * wireInterval));

    m_isNormalised = false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool GridManager::Grid::IsInsideGrid(const TVector3 &position, const float width) const
{
    const float gridMinWireCoord(std::min(m_wireBoundaries.front(), m_wireBoundaries.back()));
    const float gridMaxWireCoord(std::max(m_wireBoundaries.front(), m_wireBoundaries.back()));
    const float hitWireCoord(position.Z());

    if (((gridMaxWireCoord - hitWireCoord) < std::numeric_limits<float>::epsilon()) ||
        ((hitWireCoord - gridMinWireCoord) < std::numeric_limits<float>::epsilon()))
    {
        return false;
    }

    const float gridMinDriftCoord(std::min(m_driftBoundaries.front(), m_driftBoundaries.back()));
    const float gridMaxDriftCoord(std::max(m_driftBoundaries.front(), m_driftBoundaries.back()));
    const float hitMinDriftCoord = (width < std::numeric_limits<float>::epsilon()) ? position.X() : position.X() - (m_nSigmaConsidered * (width / 2.0));
    const float hitMaxDriftCoord = (width < std::numeric_limits<float>::epsilon()) ? position.X() : position.X() + (m_nSigmaConsidered * (width / 2.0));

    if (((gridMaxDriftCoord - hitMinDriftCoord) < std::numeric_limits<float>::epsilon()) ||
        ((hitMaxDriftCoord - gridMinDriftCoord) < std::numeric_limits<float>::epsilon()))
    {
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GridManager::Grid::AddToGrid(const TVector3 &position, const float width, const float energy)
{
    if (width < std::numeric_limits<float>::epsilon())
        return;
    
    const float wireInterval = std::fabs(m_wireBoundaries.at(0) - m_wireBoundaries.at(1));
    const int wireBin = (m_wireBoundaries.back() > m_wireBoundaries.front()) ?
        std::floor((position.Z() - m_wireBoundaries.front()) / wireInterval) :
        std::floor((m_wireBoundaries.front() - position.Z()) / wireInterval);

    if ((wireBin < 0) || (wireBin >= static_cast<int>(m_axisDimensions)))
        return;

    // Now fill assuming hits are Gaussian...
    const float driftInterval = std::fabs(m_driftBoundaries.at(0) - m_driftBoundaries.at(1));
    const float hitLowEdge = position.X() - (m_nSigmaConsidered * (width / 2.f));
    const float hitHighEdge = position.X() + (m_nSigmaConsidered * (width / 2.f));

    // Link up 'low X' and 'high X' with grid's 'start' and 'end' definitions
    float hitStartEdge = hitLowEdge;
    float hitEndEdge = hitHighEdge;
    int startDriftBin = std::floor((hitStartEdge - m_driftBoundaries.front()) / driftInterval); 
    int endDriftBin = std::floor((hitEndEdge - m_driftBoundaries.front()) / driftInterval);

    if (m_driftBoundaries.back() < m_driftBoundaries.front())
    {
        hitStartEdge = hitHighEdge;
        hitEndEdge = hitLowEdge;
        startDriftBin = std::floor((m_driftBoundaries.front() - hitStartEdge) / driftInterval);
        endDriftBin = std::floor((m_driftBoundaries.front() - hitEndEdge) / driftInterval);
    }

    // Loop over the drift bins, and fill grid
    for (int iDriftBin = startDriftBin; iDriftBin <= endDriftBin; ++iDriftBin)
    {
        if (iDriftBin < 0)
            continue;

        if (iDriftBin >= static_cast<int>(m_axisDimensions))
            continue;

        const float integralStartX = (iDriftBin == startDriftBin) ? hitStartEdge : m_driftBoundaries.at(iDriftBin);
        const float integralEndX = (iDriftBin == endDriftBin) ? hitEndEdge : m_driftBoundaries.at(iDriftBin + 1);

        // Integrate the Gaussian area...
        const float chargeFraction = IvysaurusUtils::IntegrateGaussian(integralStartX, integralEndX, position.X(), (width / 2.f), m_integralStep); 
        const float entryEnergy = energy * chargeFraction;

        // Now fill grid
        m_gridValues[iDriftBin][wireBin].first += entryEnergy;
        m_gridValues[iDriftBin][wireBin].second = true;
    }    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GridManager::Grid::NormaliseGrid(const float mean, const float std)
{
    if (m_isNormalised)
        throw cet::exception("ivysaurus::Grid") << "the entries are already normalised!";

    if (std < std::numeric_limits<float>::epsilon())
        throw cet::exception("ivysaurus::Grid") << "can't normalise with a std equal to 0!";

    for (unsigned int driftIndex = 0; driftIndex < m_axisDimensions; ++driftIndex)
    {
        for (unsigned int wireIndex = 0; wireIndex < m_axisDimensions; ++wireIndex)
        {
            // Only normalise if entry is filled...
            if (!m_gridValues[driftIndex][wireIndex].second)
                continue;

            float &gridEntry = m_gridValues[driftIndex][wireIndex].first;
            gridEntry = std::log1p(gridEntry);
            gridEntry = ((gridEntry - mean) / std);
        }
    }

    m_isNormalised = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

GridManager::GridManager(const fhicl::ParameterSet& pset) :
    m_hitModuleLabel(pset.get<std::string>("HitModuleLabel")),
    m_recoModuleLabel(pset.get<std::string>("RecoModuleLabel")),
    m_trackModuleLabel(pset.get<std::string>("TrackModuleLabel")),
    m_showerModuleLabel(pset.get<std::string>("ShowerModuleLabel")),
    m_shrEndpointLengthFrac(pset.get<float>("ShrEndpointLengthFrac")),
    m_gridSize3D(pset.get<float>("GridSize3D")),
    m_dimensions(pset.get<unsigned int>("GridDimensions")),
    m_gridMean(pset.get<float>("GridMean")),
    m_gridStd(pset.get<float>("GridStd")),
    m_nSigmaConsidered(pset.get<unsigned int>("NSigmaConsidered")),
    m_integralStep(pset.get<float>("IntegralStep")),
    m_recombFactor(pset.get<float>("RecombFactor")),
    m_calorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

GridManager::~GridManager()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

GridManager::GridMap GridManager::ObtainGridMap(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, const bool isStart) const
{
    GridManager::GridMap gridMap;

    // Get 3D grid center
    TVector3 center3D(0.f, 0.f, 0.f);
    if (!this->GetGridCenter(evt, pfparticle, isStart, center3D))
        return gridMap;

    // Now create grids        
    for (IvysaurusUtils::PandoraView pandoraView : {IvysaurusUtils::PandoraView::TPC_VIEW_U, 
         IvysaurusUtils::PandoraView::TPC_VIEW_V, IvysaurusUtils::PandoraView::TPC_VIEW_W})
    {
        const TVector3 center2D = ProjectIntoPandoraView(center3D, pandoraView);
        const float halfGrid = 0.5f * m_gridSize3D;
        const TVector3 bottomLHSCorner = TVector3(center2D.X() - halfGrid,
                                                  0.f,
                                                  center2D.Z() - halfGrid);
        
        gridMap.insert(std::make_pair(pandoraView, Grid(bottomLHSCorner, m_gridSize3D, m_gridSize3D, 
            m_dimensions, m_nSigmaConsidered, m_integralStep)));
    }

    std::vector<art::Ptr<recob::SpacePoint>> spToConsider;
    this->GetSpacePointsToConsider(evt, pfparticle, center3D, spToConsider);
    this->FillGrids(evt, spToConsider, gridMap);

    return gridMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------    

bool GridManager::GetGridCenter(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, const bool isStart,
    TVector3 &center) const
{
    if (isStart)
    {
        return GetStart(evt, pfparticle, center);
    }
    else
    {
        if (pfparticle->PdgCode() == pandora::MU_MINUS)
            return GetEndTrack(evt, pfparticle, center);
        else
            return GetEndShower(evt, pfparticle, center);
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool GridManager::GetStart(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, TVector3 &start) const
{
    if (!lar_pandora::PandoraPFParticleUtils::HasVertex(pfparticle, evt, m_recoModuleLabel))
        return false;
    
    const art::Ptr<recob::Vertex> vertex = lar_pandora::PandoraPFParticleUtils::GetVertex(pfparticle, evt, m_recoModuleLabel);
    start = TVector3(vertex->position().X(), vertex->position().Y(), vertex->position().Z());    
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool GridManager::GetEndTrack(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, TVector3 &end) const
{
    if (!lar_pandora::PandoraPFParticleUtils::HasTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel))
        return false;

    const art::Ptr<recob::Track> track = lar_pandora::PandoraPFParticleUtils::GetTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel);
    end = TVector3(track->End().X(), track->End().Y(), track->End().Z());
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool GridManager::GetEndShower(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, TVector3 &end) const
{
    if (!lar_pandora::PandoraPFParticleUtils::HasShower(pfparticle, evt, m_recoModuleLabel, m_showerModuleLabel))
        return false;

    const art::Ptr<recob::Shower> shower = lar_pandora::PandoraPFParticleUtils::GetShower(pfparticle, evt, m_recoModuleLabel, m_showerModuleLabel);
    const TVector3 direction = TVector3(shower->Direction().X(), shower->Direction().Y(), shower->Direction().Z());
    const TVector3 start = TVector3(shower->ShowerStart().X(), shower->ShowerStart().Y(), shower->ShowerStart().Z());
    const float length = (shower->Length() * m_shrEndpointLengthFrac); // ATTN: Don't put it at the VERY end... 
    end = start + (length * direction);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GridManager::GetSpacePointsToConsider(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle,
    const TVector3 &center, std::vector<art::Ptr<recob::SpacePoint>> &spToConsider) const
{
    std::vector<art::Ptr<recob::SpacePoint>> spacepoints = lar_pandora::PandoraPFParticleUtils::GetSpacePoints(pfparticle, evt, m_recoModuleLabel);
    this->GetSpacePointsToConsider(center, spacepoints, spToConsider);

    // Add in children spacepoints...
    const std::vector<art::Ptr<recob::PFParticle>> pfpChildren = lar_pandora::PandoraPFParticleUtils::GetChildParticles(pfparticle, evt, m_recoModuleLabel);
    for (const art::Ptr<recob::PFParticle> &childPFP : pfpChildren)
    {
        const std::vector<art::Ptr<recob::SpacePoint>> &childSpacepoints = lar_pandora::PandoraPFParticleUtils::GetSpacePoints(childPFP, evt, m_recoModuleLabel);
        this->GetSpacePointsToConsider(center, childSpacepoints, spToConsider);
    }

    // Add in parent
    unsigned int parentSelf(pfparticle->Parent());
    const std::vector<art::Ptr<recob::PFParticle>> &pfps = lar_pandora::PandoraEventUtils::GetPFParticles(evt, m_recoModuleLabel);
    for (art::Ptr<recob::PFParticle> pfp : pfps)
    {
        if (pfp->Self() == parentSelf)
        {
            if ((std::abs(pfp->PdgCode()) != pandora::NU_E) && (std::abs(pfp->PdgCode()) != pandora::NU_MU) && (std::abs(pfp->PdgCode()) != pandora::NU_TAU))
            {
                const std::vector<art::Ptr<recob::SpacePoint>> &parentSpacepoints = lar_pandora::PandoraPFParticleUtils::GetSpacePoints(pfp, evt, m_recoModuleLabel);
                this->GetSpacePointsToConsider(center, parentSpacepoints, spToConsider);
            }
            
            break;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GridManager::GetSpacePointsToConsider(const TVector3 &center3D, const std::vector<art::Ptr<recob::SpacePoint>> &pfpSPs,
    std::vector<art::Ptr<recob::SpacePoint>> &spToConsider) const
{
    for (const art::Ptr<recob::SpacePoint> &spacepoint : pfpSPs)
    {
        if (std::fabs(spacepoint->XYZ()[0] - center3D.X()) < (m_gridSize3D * 0.5f))
            spToConsider.emplace_back(spacepoint);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GridManager::FillGrids(const art::Event &evt, const std::vector<art::Ptr<recob::SpacePoint>> &spToConsider, GridManager::GridMap &gridMap) const
{
    for (const art::Ptr<recob::SpacePoint> &spacepoint : spToConsider)
    {
        const std::vector<art::Ptr<recob::Hit>> hits = lar_pandora::PandoraSpacePointUtils::GetHits(spacepoint, evt, m_recoModuleLabel);
        if (hits.empty()) { continue; }
        const art::Ptr<recob::Hit> hit = hits.front();
        const IvysaurusUtils::PandoraView thisPandoraView = IvysaurusUtils::GetPandora2DView(hit);
        GridManager::Grid &grid = gridMap.at(thisPandoraView);
        
        float hitWidth = 0.f;
        TVector3 pandoraHitPosition = TVector3(0.f, 0.f, 0.f);
        IvysaurusUtils::ObtainPandoraHitPositionAndWidth(evt, hit, thisPandoraView, pandoraHitPosition, hitWidth);

        if (!grid.IsInsideGrid(pandoraHitPosition, hitWidth))
            continue;

        const float energy = ObtainHitEnergy(evt, hit);
        grid.AddToGrid(pandoraHitPosition, hitWidth, energy);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float GridManager::ObtainHitEnergy(const art::Event &evt, const art::Ptr<recob::Hit> &hit) const
{
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

    const double charge = lar_pandora::PandoraHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, {hit});
    const double nElectrons = m_calorimetryAlg.ElectronsFromADCArea(charge, hit->WireID().Plane);
    const double hitEnergy = nElectrons / m_recombFactor / util::kGeVToElectrons;

    return hitEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GridManager::NormaliseGrid(GridManager::Grid &grid)
{
    grid.NormaliseGrid(m_gridMean, m_gridStd);
}

} //namespace ivysaurus

