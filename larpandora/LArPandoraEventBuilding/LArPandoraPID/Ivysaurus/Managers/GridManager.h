////////////////////////////////////////////////////////////////////////
/// \file    GridManager.cxx
/// \brief   A class to manage the Ivysaurus 2D grid input 
/// \author  Isobel Mawby - i.mawby1@lancaster.ac.uk
////////////////////////////////////////////////////////////////////////

#ifndef GRIDMANAGER_H
#define GRIDMANAGER_H

#include <vector>
#include <string>
#include <map>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"

#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Utils/IvysaurusUtils.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"

namespace ivysaurus
{

class GridManager
{
  public:
  class Grid
  {
    public:
      Grid(const TVector3 origin, const float driftSpan, const float wireSpan, const unsigned int dimensions, 
          const unsigned int nSigmaConsidered, const float integralStep);

      unsigned int GetAxisDimensions() const;
      std::vector<float> GetDriftBoundaries() const;
      std::vector<float> GetWireBoundaries() const;
      //float GetGridEntry(const unsigned int driftIndex, const unsigned int wireIndex) const;
      //void SetGridEntry(const unsigned int driftIndex, const unsigned int wireIndex, const float value);
      std::vector<std::vector<std::pair<float, bool>>> GetGridValues() const;
      bool IsNormalised() const;
      bool IsInsideGrid(const TVector3 &position, const float width) const;
      void AddToGrid(const TVector3 &position, const float width, const float energy);
      void NormaliseGrid(const float mean, const float std);

    private:
      unsigned int m_axisDimensions; // number of bins on each axis
      unsigned int m_nSigmaConsidered;
      float m_integralStep; 
      std::vector<float> m_driftBoundaries;
      std::vector<float> m_wireBoundaries;
      std::vector<std::vector<std::pair<float, bool>>> m_gridValues; // driftBin : [wireBins]
      bool m_isNormalised;
  };

    typedef std::map<IvysaurusUtils::PandoraView, Grid> GridMap;

    GridManager(const fhicl::ParameterSet& pset);
    ~GridManager();

    // Function to place the grid in space
    GridMap ObtainGridMap(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        const bool isStart) const;

    void NormaliseGrid(GridManager::Grid &grid);

  private:
    bool GetGridExtremalPoints(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, const bool isStart,
        TVector3 &position1, TVector3 &position2) const;
    
    bool GetStartExtremalPoints(const art::Event &evt, const std::vector<art::Ptr<recob::SpacePoint>> &spacepointsToConsider,
        const art::Ptr<recob::PFParticle> &pfparticle, TVector3 &position1, TVector3 &position2) const;

    bool GetEndExtremalPointsTrack(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        TVector3 &position1, TVector3 &position2) const;

    bool GetEndExtremalPointsShower(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        TVector3 &position1, TVector3 &position2) const;

    void GetSpacePointsToConsider(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle,
        const TVector3 &position1, const TVector3 &position2, std::vector<art::Ptr<recob::SpacePoint>> &spToConsider) const;

    void FillGrids(const art::Event &evt, const std::vector<art::Ptr<recob::SpacePoint>> &spToConsider, GridManager::GridMap &gridMap) const;
    
    float ObtainHitEnergy(const art::Event &evt, const art::Ptr<recob::Hit> &hit) const;

    std::string m_hitModuleLabel;
    std::string m_recoModuleLabel;
    std::string m_trackModuleLabel;
    std::string m_showerModuleLabel;

    float m_gridSize3D;
    unsigned int m_dimensions;
    float m_gridMean;
    float m_gridStd;
    unsigned int m_nSigmaConsidered;                                                                                                                                                                                           
    float m_integralStep;
    float m_recombFactor;
    calo::CalorimetryAlg m_calorimetryAlg;
};

/////////////////////////////////////////////////////////////

inline unsigned int GridManager::Grid::GetAxisDimensions() const 
{ 
    return m_axisDimensions; 
}

/////////////////////////////////////////////////////////////

inline std::vector<float> GridManager::Grid::GetDriftBoundaries() const 
{ 
    return m_driftBoundaries; 
}

/////////////////////////////////////////////////////////////

inline std::vector<float> GridManager::Grid::GetWireBoundaries() const 
{ 
    return m_wireBoundaries; 
}

/////////////////////////////////////////////////////////////

// inline float GridManager::Grid::GetGridEntry(const unsigned int driftIndex, 
//     const unsigned int wireIndex) const 
// { 
//     return m_gridValues[driftIndex][wireIndex].first; 
// }

/////////////////////////////////////////////////////////////

// inline void GridManager::Grid::SetGridEntry(const unsigned int driftIndex, const unsigned int wireIndex, 
//     const float value)
// { 
//     m_gridValues[driftIndex][wireIndex].first = value;
//     m_gridValues[driftIndex][wireIndex].second = true;    
// }

/////////////////////////////////////////////////////////////

inline std::vector<std::vector<std::pair<float, bool>>> GridManager::Grid::GetGridValues() const
{ 
    return m_gridValues; 
}

/////////////////////////////////////////////////////////////

inline bool GridManager::Grid::IsNormalised() const 
{ 
    return m_isNormalised;
}

}

#endif  // GRIDMANAGER_H
