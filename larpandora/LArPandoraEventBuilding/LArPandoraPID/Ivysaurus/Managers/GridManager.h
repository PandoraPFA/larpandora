/**
 *  @file  larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/GridManager.h
 *
 *  @brief A class to manage the Ivysaurus 2D grid input 
 *
 */
#ifndef GRIDMANAGER_H
#define GRIDMANAGER_H
// ART
#include "art/Framework/Principal/Event.h"
// LArSoft
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Utils/IvysaurusUtils.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

namespace ivysaurus
{

/**
 *  @brief  GridManager class
 */    
class GridManager
{
  public:
      typedef std::vector<std::vector<std::pair<float, bool>>> GridValues;
    
      /**
       *  @brief  Grid class
       */    
      class Grid
      {
      public:
         /**
          *  @brief  Default constructor
          */      
         Grid(const TVector3 origin, const float driftSpan, const float wireSpan, const unsigned int dimensions, 
              const unsigned int nSigmaConsidered, const float integralStep);

          /**
           *  @brief  Get the dimensions of the square grid
           *
           *  @return the grid dimensions
           */
          unsigned int GetAxisDimensions() const;

          /**
           *  @brief  Get the bin boundaries for the drift axis
           *
           *  @return the drift axis bin boundaries
           */          
          std::vector<float> GetDriftBoundaries() const;

          /**
           *  @brief  Get the bin boundaries for the wire axis
           *
           *  @return the wire axis bin boundaries
           */                    
          std::vector<float> GetWireBoundaries() const;

          /**
           *  @brief  Get the grid entries
           *
           *  @return the grid entries
           */          
          GridManager::GridValues GetGridValues() const;

          /**
           *  @brief  Return whether the grid has been normalised
           *
           *  @return whether the grid has been normalised
           */          
          bool IsNormalised() const;

          /**
           *  @brief  Determine whether an input position is inside the grid (function consideres hit width)
           *
           *  @param position the input hit position
           *  @param width the input hit width
           *
           *  @return whether the hit lies within the grid
           */                    
          bool IsInsideGrid(const TVector3 &position, const float width) const;

           /**
            *  @brief  Add an input hit's energy to the grid
            *
            *  @param position the input hit position
            *  @param width the input hit width
            *  @param energy the input hit energy (corrected)
            */
          void AddToGrid(const TVector3 &position, const float width, const float energy);

           /**
            *  @brief  Normalise the grid entries
            *
            *  @param mean the mean of the grid entries (calculated offline)
            *  @param std the standard deviation of the grid entries (calculated offline)
            */
          void NormaliseGrid(const float mean, const float std);

      private:
          unsigned int m_axisDimensions;         ///< the dimensions of the square grid
          unsigned int m_nSigmaConsidered;       ///< the sigma range of hit energy considered when grid filling
          float m_integralStep;                  ///< the step size for numerical integration of hit energy
          std::vector<float> m_driftBoundaries;  ///< the bin boundaries for the drift axis
          std::vector<float> m_wireBoundaries;   ///< the bin boundaries for the wire axis
          GridManager::GridValues m_gridValues;  ///< the grid entries driftBin->wireBin->[energy, isFilled]
          bool m_isNormalised;                   ///< whether the grid entried have are normalised
  };

    typedef std::map<IvysaurusUtils::PandoraView, Grid> GridMap;
    
    /**
     *  @brief  Default constructor
     */
    GridManager(const fhicl::ParameterSet& pset);

    /**
     *  @brief  Default destructor
     */    
    ~GridManager();

    /**
     *  @brief  Obtain the grid map (view->grid) for a given PFParticle
     *          ATTN: returned grids are not normalised
     *
     *  @param evt the art event
     *  @param pfparticle the input PFParticle
     *  @param isStart whether to calculate the grids at the start or end of the PFParticle
     *
     *  @return the grid map for the input PFParticle
     */
    GridMap ObtainGridMap(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        const bool isStart) const;

    /**
     *  @brief  Normalise the FILLED grid entries, unfilled entries are not normalised
     *
     *  @param grid the input grid
     */
    void NormaliseGrid(GridManager::Grid &grid);

  private:
    /**
     *  @brief  Use track/shower parameterisations to obtain the boundaries of the region,
     *          along particle's 3D trajectory, that will be used to fill the 2D grids
     *
     *  @param evt the art event
     *  @param pfparticle the input PFParticle
     *  @param isStart whether to calculate the grids at the start or end of the PFParticle
     *  @param[out] center the grid center (in 3D)
     *
     *  @return whether the extremal points could be calculated (sometimes fits/direction estimators fail)
     */
    bool GetGridCenter(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, const bool isStart,
        TVector3 &center) const;

    /**
     *  @brief  Set the central grid position as the PFParticle vertex
     *
     *  @param evt the art event
     *  @param pfparticle the input PFParticle
     *  @param[out] center the grid center (in 3D)
     *
     *  @return whether the extremal points could be calculated (sometimes vertexing fails)
     */
    bool GetStart(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, TVector3 &center) const;

    /**
     *  @brief  Set the central grid position as the Track endpoint
     *
     *  @param evt the art event
     *  @param pfparticle the input PFParticle
     *  @param[out] end the grid center (in 3D)
     *
     *  @return whether the endpoint could be calculated (sometimes track fits fail)
     */
    bool GetEndTrack(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, TVector3 &end) const;

    /**
     *  @brief  Set the central grid position as the shifted Track endpoint       
     *
     *  @param evt the art event
     *  @param pfparticle the input PFParticle
     *  @param[out] end the grid center (in 3D)
     *
     *  @return whether the extremal points could be calculated (sometimes track fits fail)
     */
    bool GetEndShower(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, TVector3 &end) const;

    /**
     *  @brief  Obtain the spacepoints associated with the PFParticle (incl. parent, and children) that
     *          potentially live in the grid (apply some generous cuts)
     *
     *  @param evt the art event
     *  @param pfparticle the input PFParticle
     *  @param center the grid center (in 3D)
     *  @param[out] spToConsider the spacepoints to consider for grid filling
     */
    void GetSpacePointsToConsider(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle,
        const TVector3 &center3D, std::vector<art::Ptr<recob::SpacePoint>> &spToConsider) const;

    /**
     *  @brief  Identify the spacepoints, from an input list, that potentially live in the grid (apply some generous cuts)
     *
     *  @param center the grid center (in 3D)
     *  @param pfpSPs input list of spacepoints to filter
     *  @param[out] spToConsider the spacepoints to consider for grid filling
     */
    void GetSpacePointsToConsider(const TVector3 &center3D, const std::vector<art::Ptr<recob::SpacePoint>> &pfpSPs,
        std::vector<art::Ptr<recob::SpacePoint>> &spToConsider) const;
    
    /**
     *  @brief  Fill the grids with the energy of the hits associated with the input spacepoints
     *
     *  @param evt the art event
     *  @param spToConsider the spacepoints to consider for grid filling
     *  @param[out] gridMap the grid map to fill
     */
    void FillGrids(const art::Event &evt, const std::vector<art::Ptr<recob::SpacePoint>> &spToConsider, GridManager::GridMap &gridMap) const;

    /**
     *  @brief  Obtain the energy of a hit, corrected for recombination effects
     *
     *  @param evt the art event
     *  @param hit the input hit
     *
     *  @return the corrected hit energy
     */
    float ObtainHitEnergy(const art::Event &evt, const art::Ptr<recob::Hit> &hit) const;

    std::string m_hitModuleLabel;    ///< the label for the hit producer
    std::string m_recoModuleLabel;   ///< the label of pandora
    std::string m_trackModuleLabel;  ///< the label for the track producer
    std::string m_showerModuleLabel; ///< the label for the shower producer
    float m_shrEndpointLengthFrac;   ///< the shower length fraction that defines the shower endpoint
    float m_gridSize3D;              ///< the 3D length of the region considered for grid filling (the grids are also square with this length)
    unsigned int m_dimensions;       ///< the number of bins along each axis of the square grid
    float m_gridMean;                ///< the mean used to normalise grid entries (calculated offline)
    float m_gridStd;                 ///< the standard deviation used to normalise grid entries (calculated offline)
    unsigned int m_nSigmaConsidered; ///< the number of sigma of hit energy considered when filling the grid
    float m_integralStep;            ///< the step size for numerical integration of hit energy
    float m_recombFactor;            ///< the recombination factor used to correct hit energy
    calo::CalorimetryAlg m_calorimetryAlg; ///< the calorimetry algorithm used to correct hit energy
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int GridManager::Grid::GetAxisDimensions() const 
{ 
    return m_axisDimensions; 
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::vector<float> GridManager::Grid::GetDriftBoundaries() const 
{ 
    return m_driftBoundaries; 
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::vector<float> GridManager::Grid::GetWireBoundaries() const 
{ 
    return m_wireBoundaries; 
}

//------------------------------------------------------------------------------------------------------------------------------------------    

inline GridManager::GridValues GridManager::Grid::GetGridValues() const
{ 
    return m_gridValues; 
}

//------------------------------------------------------------------------------------------------------------------------------------------    

inline bool GridManager::Grid::IsNormalised() const 
{ 
    return m_isNormalised;
}
    
} // namepace ivysaurus

#endif  // GRIDMANAGER_H

