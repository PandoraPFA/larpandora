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
      Grid(const TVector3 origin, const float driftSpan, const float wireSpan, const unsigned int dimensions, const float maxGridEntry, 
          const unsigned int nSigmaConsidered, const float integralStep, const IvysaurusUtils::PandoraView pandoraView, const bool isInitialised);

      unsigned int GetAxisDimensions() const;
      std::vector<float> GetDriftBoundaries() const;
      std::vector<float> GetWireBoundaries() const;
      float GetGridEntry(const unsigned int driftIndex, const unsigned int wireIndex) const;
      void SetGridEntry(const unsigned int driftIndex, const unsigned int wireIndex, const float value);
      void ResetGridEntries();
      std::vector<std::vector<float>> GetGridValues() const;
      IvysaurusUtils::PandoraView GetPandoraView() const;
      bool IsInitialised() const;
      bool IsNormalised() const;
      bool IsInsideGrid(const TVector3 &position, const float width) const;
      void AddToGrid(const TVector3 &position, const float width, const float energy);
      void NormaliseGrid();

    private:
      unsigned int m_axisDimensions; // number of bins on each axis
      float m_maxGridEntry;
      unsigned int m_nSigmaConsidered;
      float m_integralStep; 
      std::vector<float> m_driftBoundaries;
      std::vector<float> m_wireBoundaries;
      std::vector<art::Ptr<recob::Hit>> m_gridHitList;
      std::vector<std::vector<float>> m_gridValues; // driftBin : [wireBins]
      IvysaurusUtils::PandoraView m_pandoraView;
      bool m_isInitialised;
      bool m_isAveraged;
      bool m_isNormalised;
  };

    typedef std::map<IvysaurusUtils::PandoraView, Grid> GridMap;

    GridManager(const fhicl::ParameterSet& pset);
    ~GridManager();

    // Function to place the grid in space
    GridMap ObtainGridMap(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        const bool isStart) const;

    void FillGrids(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, GridManager::GridMap &gridMap) const;    

  private:
    bool GetStartExtremalPoints(const art::Event &evt, const std::vector<art::Ptr<recob::SpacePoint>> &spacepointsToConsider,
        const art::Ptr<recob::PFParticle> &pfparticle, TVector3 &position1, TVector3 &position2) const;

    bool GetEndExtremalPointsTrack(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        TVector3 &position1, TVector3 &position2) const;

    bool GetEndExtremalPointsShower(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        TVector3 &position1, TVector3 &position2) const;

    float ObtainHitEnergy(const art::Event &evt, const art::Ptr<recob::Hit> &hit) const;

    std::string m_hitModuleLabel;
    std::string m_recoModuleLabel;
    std::string m_trackModuleLabel;
    std::string m_showerModuleLabel;

    float m_gridSize3D;
    unsigned int m_dimensions;
    bool m_addChildrenToGrid;
    float m_maxGridEntry;
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

inline float GridManager::Grid::GetGridEntry(const unsigned int driftIndex, 
    const unsigned int wireIndex) const 
{ 
    return m_gridValues[driftIndex][wireIndex]; 
}

/////////////////////////////////////////////////////////////

inline void GridManager::Grid::SetGridEntry(const unsigned int driftIndex, const unsigned int wireIndex, 
    const float value)
{ 
    m_gridValues[driftIndex][wireIndex] = value;
}

/////////////////////////////////////////////////////////////

inline void GridManager::Grid::ResetGridEntries()
{ 
    m_gridValues = std::vector<std::vector<float>>(m_axisDimensions, std::vector<float>(m_axisDimensions, 0.0));
    m_isNormalised = false;
}

/////////////////////////////////////////////////////////////

inline std::vector<std::vector<float>> GridManager::Grid::GetGridValues() const 
{ 
    return m_gridValues; 
}

/////////////////////////////////////////////////////////////

inline IvysaurusUtils::PandoraView GridManager::Grid::GetPandoraView() const 
{ 
    return m_pandoraView; 
}

/////////////////////////////////////////////////////////////

inline bool GridManager::Grid::IsInitialised() const 
{ 
    return m_isInitialised;
}

/////////////////////////////////////////////////////////////

inline bool GridManager::Grid::IsNormalised() const 
{ 
    return m_isNormalised;
}

}

#endif  // GRIDMANAGER_H
