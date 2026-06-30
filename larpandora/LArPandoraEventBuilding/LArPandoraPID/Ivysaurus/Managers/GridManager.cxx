////////////////////////////////////////////////////////////////////////
/// \file    GridManager.cxx
/// \brief   A class to manage the Ivysaurus 2D grid input 
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

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"

#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"

#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/GridManager.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Utils/IvysaurusUtils.h"
#include "larpandora/LArPandoraUtils/PandoraEventUtils.h"
#include "larpandora/LArPandoraUtils/PandoraPFParticleUtils.h"
#include "larpandora/LArPandoraUtils/PandoraSpacePointUtils.h"
#include "larpandora/LArPandoraUtils/PandoraHitUtils.h"


namespace ivysaurus
{

GridManager::Grid::Grid(const TVector3 origin, const float driftSpan, const float wireSpan, 
    const unsigned int dimensions,const unsigned int nSigmaConsidered, const float integralStep) : 
        m_axisDimensions(dimensions),
        m_nSigmaConsidered(nSigmaConsidered),
        m_integralStep(integralStep)
{
    m_gridValues = std::vector<std::vector<std::pair<float, bool>>>(m_axisDimensions, std::vector<std::pair<float, bool>>(m_axisDimensions, std::make_pair(0.f, false)));    

    const float driftInterval = driftSpan / m_axisDimensions;

    for (unsigned int i = 0; i <= m_axisDimensions; ++i)
        m_driftBoundaries.push_back(origin.X() + (i * driftInterval));

    const float wireInterval = wireSpan / m_axisDimensions;

    for (unsigned int i = 0; i <= m_axisDimensions; ++i)
        m_wireBoundaries.push_back(origin.Z() + (i * wireInterval));

    m_isNormalised = false;
}

/////////////////////////////////////////////////////////////

bool GridManager::Grid::IsInsideGrid(const TVector3 &position, const float width) const
{
    // Wire axis
    const float gridMinWireCoord(std::min(m_wireBoundaries.front(), m_wireBoundaries.back()));
    const float gridMaxWireCoord(std::max(m_wireBoundaries.front(), m_wireBoundaries.back()));
    const float hitWireCoord(position.Z());

    if (std::fabs(gridMinWireCoord - hitWireCoord) < std::numeric_limits<float>::epsilon())
        return false;

    if (std::fabs(gridMaxWireCoord - hitWireCoord) < std::numeric_limits<float>::epsilon())
        return false;

    if ((hitWireCoord < gridMinWireCoord) || (hitWireCoord > gridMaxWireCoord))
        return false;

    // Drift axis
    const float gridMinDriftCoord(std::min(m_driftBoundaries.front(), m_driftBoundaries.back()));
    const float gridMaxDriftCoord(std::max(m_driftBoundaries.front(), m_driftBoundaries.back()));
    const float hitMinDriftCoord = (width < std::numeric_limits<float>::epsilon()) ? position.X() : position.X() - (m_nSigmaConsidered * (width / 2.0));
    const float hitMaxDriftCoord = (width < std::numeric_limits<float>::epsilon()) ? position.X() : position.X() + (m_nSigmaConsidered * (width / 2.0));

    if (std::fabs(hitMaxDriftCoord - gridMinDriftCoord) < std::numeric_limits<float>::epsilon())
        return false;

    if (std::fabs(hitMinDriftCoord - gridMaxDriftCoord) < std::numeric_limits<float>::epsilon())
        return false;

    if ((hitMaxDriftCoord < gridMinDriftCoord) || (hitMinDriftCoord > gridMaxDriftCoord))
        return false;

    return true;
}

/////////////////////////////////////////////////////////////

void GridManager::Grid::AddToGrid(const TVector3 &position, const float width, const float energy)
{
    // Get wire bin
    const float wireInterval = std::fabs(m_wireBoundaries.at(0) - m_wireBoundaries.at(1));
    int wireBin = std::floor((position.Z() - m_wireBoundaries.front()) / wireInterval); 

    if (m_wireBoundaries.back() < m_wireBoundaries.front())
        wireBin = std::floor((m_wireBoundaries.front() - position.Z()) / wireInterval);

    if (wireBin < 0)
        return;

    if (wireBin >= static_cast<int>(m_axisDimensions))
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

/////////////////////////////////////////////////////////////

void GridManager::Grid::NormaliseGrid(const float mean, const float std)
{
    if (m_isNormalised)
        throw cet::exception("ivysaur::Grid") << "the entries are already normalised!";

    for (unsigned int driftIndex = 0; driftIndex < m_axisDimensions; ++driftIndex)
    {
        for (unsigned int wireIndex = 0; wireIndex < m_axisDimensions; ++wireIndex)
        {
            // Only normalise if entry is filled...
            if (!m_gridValues[driftIndex][wireIndex].second)
                continue;

            float gridEntry = m_gridValues[driftIndex][wireIndex].first;
            gridEntry = std::log1p(gridEntry);
            m_gridValues[driftIndex][wireIndex].first = ((gridEntry - mean) / std);
        }
    }

    m_isNormalised = true;
}


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

GridManager::GridManager(const fhicl::ParameterSet& pset) :
    m_hitModuleLabel(pset.get<std::string>("HitModuleLabel")),
    m_recoModuleLabel(pset.get<std::string>("RecoModuleLabel")),
    m_trackModuleLabel(pset.get<std::string>("TrackModuleLabel")),
    m_showerModuleLabel(pset.get<std::string>("ShowerModuleLabel")),
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

/////////////////////////////////////////////////////////////

GridManager::~GridManager()
{
}


/////////////////////////////////////////////////////////////

GridManager::GridMap GridManager::ObtainGridMap(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, const bool isStart) const
{
    GridManager::GridMap gridMap;

    // Get 3D bounding box
    TVector3 position1(0.f, 0.f, 0.f), position2(0.f, 0.f, 0.f);    
    if (!this->GetGridExtremalPoints(evt, pfparticle, isStart, position1, position2))
        return gridMap;

    // Collect spacepoints
    std::vector<art::Ptr<recob::SpacePoint>> spToConsider;
    this->GetSpacePointsToConsider(evt, pfparticle, position1, position2, spToConsider);

    // Now create grids
    for (IvysaurusUtils::PandoraView pandoraView : {IvysaurusUtils::PandoraView::TPC_VIEW_U, 
         IvysaurusUtils::PandoraView::TPC_VIEW_V, IvysaurusUtils::PandoraView::TPC_VIEW_W})
    {
        const TVector3 projectedPosition1 = ProjectIntoPandoraView(position1, pandoraView);
        const TVector3 projectedSeed = ProjectIntoPandoraView(position2, pandoraView);

        float driftScale(projectedSeed.X() > projectedPosition1.X() ? 1.f : -1.f);
        float wireScale(projectedSeed.Z() > projectedPosition1.Z() ? 1.f : -1.f);

        
        const TVector3 projectedPosition2 = TVector3(projectedPosition1.X() + (m_gridSize3D * driftScale),
                                                     0.f,
                                                     projectedPosition1.Z() + (m_gridSize3D * wireScale));
        
        const float driftSpan = projectedPosition2.X() - projectedPosition1.X();
        const float wireSpan = projectedPosition2.Z() - projectedPosition1.Z();

        gridMap.insert(std::make_pair(pandoraView, Grid(projectedPosition1, driftSpan, wireSpan, 
            m_dimensions, m_nSigmaConsidered, m_integralStep)));
    }

    // Now fill hits
    this->FillGrids(evt, spToConsider, gridMap);

    return gridMap;
}
    
/////////////////////////////////////////////////////////////

bool GridManager::GetGridExtremalPoints(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, const bool isStart,
    TVector3 &position1, TVector3 &position2) const
{
    // Need spacepoints to calculate initial particle direction
    const std::vector<art::Ptr<recob::SpacePoint>> spacepoints = lar_pandora::PandoraPFParticleUtils::GetSpacePoints(pfparticle, evt, m_recoModuleLabel); 
    if (spacepoints.empty()) { return false; }

    if (isStart)
    {
        if (!GetStartExtremalPoints(evt, spacepoints, pfparticle, position1, position2))
            return false;
    }
    else
    {
        if (pfparticle->PdgCode() == 13)
        {
            if (!GetEndExtremalPointsTrack(evt, pfparticle, position1, position2))
                if (!GetEndExtremalPointsShower(evt, pfparticle, position1, position2))
                    return false;
        }
        else
        {
            if (!GetEndExtremalPointsShower(evt, pfparticle, position1, position2))
                if (!GetEndExtremalPointsTrack(evt, pfparticle, position1, position2))
                    return false;
        }
    }

    return true;
}
/////////////////////////////////////////////////////////////

bool GridManager::GetStartExtremalPoints(const art::Event &evt, const std::vector<art::Ptr<recob::SpacePoint>> &spacepointsToConsider, 
    const art::Ptr<recob::PFParticle> &pfparticle, TVector3 &position1, TVector3 &position2) const
{
    if (!lar_pandora::PandoraPFParticleUtils::HasVertex(pfparticle, evt, m_recoModuleLabel)) { return false; }
    const art::Ptr<recob::Vertex> vertex = lar_pandora::PandoraPFParticleUtils::GetVertex(pfparticle, evt, m_recoModuleLabel);
    const TVector3 vertexPos = TVector3(vertex->position().X(), vertex->position().Y(), vertex->position().Z());

    TVector3 initialDir3D(0.f, 0.f, 0.f);
    if (!IvysaurusUtils::GetInitialDirection(evt, vertexPos, spacepointsToConsider, m_recoModuleLabel, initialDir3D))
        return false;

    position1 = TVector3(vertex->position().X(), vertex->position().Y(), vertex->position().Z());
    const float diagonalLength = sqrt(2.0 * (m_gridSize3D * m_gridSize3D));
    //position2 = position1 + (initialDir3D * diagonalLength);

    position1 = position1 - (initialDir3D * diagonalLength * 0.10);
    position2 = position1 + (initialDir3D * diagonalLength * 0.90);
    
    return true;
}

/////////////////////////////////////////////////////////////

bool GridManager::GetEndExtremalPointsTrack(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    TVector3 &position1, TVector3 &position2) const
{
    // Has associated track object?
    if (!lar_pandora::PandoraPFParticleUtils::HasTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel))
        return false;

    const art::Ptr<recob::Track> track = lar_pandora::PandoraPFParticleUtils::GetTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel);

    const TVector3 direction = TVector3(track->EndDirection().X(), track->EndDirection().Y(), track->EndDirection().Z());
    const TVector3 end = TVector3(track->End().X(), track->End().Y(), track->End().Z());

    const float halfLength = 0.5f * std::sqrt(2.f) * m_gridSize3D;
    position1 = end - (direction * halfLength); 
    position2 = end + (direction * halfLength);

    return true;
}

/////////////////////////////////////////////////////////////

bool GridManager::GetEndExtremalPointsShower(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    TVector3 &position1, TVector3 &position2) const
{
    // Has associated shower object?
    if (!lar_pandora::PandoraPFParticleUtils::HasShower(pfparticle, evt, m_recoModuleLabel, m_showerModuleLabel))
        return false;

    const art::Ptr<recob::Shower> shower = lar_pandora::PandoraPFParticleUtils::GetShower(pfparticle, evt, m_recoModuleLabel, m_showerModuleLabel);

    const TVector3 direction = TVector3(shower->Direction().X(), shower->Direction().Y(), shower->Direction().Z());
    const TVector3 start = TVector3(shower->ShowerStart().X(), shower->ShowerStart().Y(), shower->ShowerStart().Z());
    const float length = (shower->Length() * 0.9f); // ATTN: Don't put it at the VERY end... 
    const TVector3 end = start + (length * direction);

    const float diagonalLength = sqrt(2.0 * (m_gridSize3D * m_gridSize3D));
    position1 = end - (direction * (diagonalLength / 2.0)); 
    position2 = end + (direction * (diagonalLength / 2.0));

    return true;
}

/////////////////////////////////////////////////////////////

void GridManager::GetSpacePointsToConsider(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle,
    const TVector3 &position1, const TVector3 &position2, std::vector<art::Ptr<recob::SpacePoint>> &spToConsider) const
{
    std::vector<art::Ptr<recob::SpacePoint>> spacepoints = lar_pandora::PandoraPFParticleUtils::GetSpacePoints(pfparticle, evt, m_recoModuleLabel);
    
    // Add in children spacepoints...
    const std::vector<art::Ptr<recob::PFParticle>> pfpChildren = lar_pandora::PandoraPFParticleUtils::GetChildParticles(pfparticle, evt, m_recoModuleLabel);
    for (const art::Ptr<recob::PFParticle> &childPFP : pfpChildren)
    {
        const std::vector<art::Ptr<recob::SpacePoint>> &childSpacepoints = lar_pandora::PandoraPFParticleUtils::GetSpacePoints(childPFP, evt, m_recoModuleLabel);
        spacepoints.insert(spacepoints.end(), childSpacepoints.begin(), childSpacepoints.end());
    }

    // Pick up spacepoints between position 1&2
    TVector3 direction = (position2 - position1).Unit();
    float mag = (position2 - position1).Mag();
    for (const art::Ptr<recob::SpacePoint> &spacepoint : spacepoints)
    {
        TVector3 displacement = spacepoint->XYZ() - position1;
        float l = direction.Dot(displacement);

        if ((l < 0.f) || (l > mag))
            continue;
    
        spToConsider.push_back(spacepoint);
    }
}

/////////////////////////////////////////////////////////////

void GridManager::FillGrids(const art::Event &evt, const std::vector<art::Ptr<recob::SpacePoint>> &spToConsider, GridManager::GridMap &gridMap) const
{
    // Add hits to grids
    for (const art::Ptr<recob::SpacePoint> &spacepoint : spToConsider)
    {
        std::vector<art::Ptr<recob::Hit>> hits = lar_pandora::PandoraSpacePointUtils::GetHits(spacepoint, evt, m_recoModuleLabel);
        if (hits.empty()) { continue; }
        art::Ptr<recob::Hit> hit = hits.front();
        const IvysaurusUtils::PandoraView thisPandoraView = IvysaurusUtils::GetPandora2DView(hit);
        GridManager::Grid &grid = gridMap.at(thisPandoraView);
        
        // Get 2D hit position and width
        float hitWidth = 0.f;
        TVector3 pandoraHitPosition = TVector3(0.f, 0.f, 0.f);
        IvysaurusUtils::ObtainPandoraHitPositionAndWidth(evt, hit, thisPandoraView, pandoraHitPosition, hitWidth);

        // Check hit is inside grid
        if (!grid.IsInsideGrid(pandoraHitPosition, hitWidth))
            continue;

        // Add its energy to the grid
        const float energy = ObtainHitEnergy(evt, hit);
        grid.AddToGrid(pandoraHitPosition, hitWidth, energy);
    }
}

/////////////////////////////////////////////////////////////

float GridManager::ObtainHitEnergy(const art::Event &evt, const art::Ptr<recob::Hit> &hit) const
{
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

    const double charge = lar_pandora::PandoraHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, {hit});
    const double nElectrons = m_calorimetryAlg.ElectronsFromADCArea(charge, hit->WireID().Plane);
    const double hitEnergy = nElectrons / m_recombFactor / util::kGeVToElectrons;

    return hitEnergy;
}

/////////////////////////////////////////////////////////////

void GridManager::NormaliseGrid(GridManager::Grid &grid)
{
    grid.NormaliseGrid(m_gridMean, m_gridStd);
}

}
