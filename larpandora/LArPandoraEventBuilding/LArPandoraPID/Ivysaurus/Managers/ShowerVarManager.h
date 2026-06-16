////////////////////////////////////////////////////////////////////////
/// \file    ShowerVarManager.h
/// \brief   A class to manage the Ivysaurus shower variables input 
/// \author  Isobel Mawby - i.mawby1@lancaster.ac.uk
////////////////////////////////////////////////////////////////////////

#ifndef SHOWERMANAGER_H
#define SHOWERMANAGER_H

#include <vector>
#include <string>
#include <map>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Utils/IvysaurusUtils.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"

namespace ivysaurus
{

class ShowerVarManager
{
  public:
  class ShowerVars
  {
    public:
      ShowerVars();

      bool GetIsNormalised() const;
      float GetDisplacement() const;
      float GetDCA() const;
      float GetTrackStubLength() const;
      float GetNuVertexAvSeparation() const;
      float GetNuVertexChargeAsymmetry() const;
      float GetFoundConnectionPathway() const;
      float GetInitialGapSize() const;
      float GetLargestGapSize() const;
      float GetPathwayLength() const;
      float GetPathwayScatteringAngle2D() const;
      float GetNShowerHits() const;
      float GetFoundHitRatio() const;
      float GetScatterAngle() const;
      float GetOpeningAngle() const;
      float GetNuVertexEnergyAsymmetry() const;
      float GetNuVertexEnergyWeightedMeanRadialDistance() const;
      float GetShowerStartEnergyAsymmetry() const;
      float GetShowerStartMoliereRadius() const;
      float GetNAmbiguousViews() const;
      float GetUnaccountedEnergy() const;
      ////
      void SetIsNormalised(const bool isNormalised);
      void SetDisplacement(const float displacement);
      void SetDCA(const float dca);
      void SetTrackStubLength(const float trackStubLength);
      void SetNuVertexAvSeparation(const float nuVertexAvSeparation);
      void SetNuVertexChargeAsymmetry(const float nuVertexChargeAsymmetry);
      void SetFoundConnectionPathway(const float foundConnectionPathway);
      void SetInitialGapSize(const float initialGapSize);
      void SetLargestGapSize(const float largestGapSize);
      void SetPathwayLength(const float pathwayLength);
      void SetPathwayScatteringAngle2D(const float pathwayScatteringAngle2D);
      void SetNShowerHits(const float nShowerHits);
      void SetFoundHitRatio(const float foundHitRatio);
      void SetScatterAngle(const float scatterAngle);
      void SetOpeningAngle(const float openingAngle);
      void SetNuVertexEnergyAsymmetry(const float nuVertexEnergyAsymmetry);
      void SetNuVertexEnergyWeightedMeanRadialDistance(const float nuVertexEnergyWeightedMeanRadialDistance);
      void SetShowerStartEnergyAsymmetry(const float showerStartEnergyAsymmetry);
      void SetShowerStartMoliereRadius(const float showerStartMoliereRadius);
      void SetNAmbiguousViews(const float nAmbiguousViews);
      void SetUnaccountedEnergy(const float unaccountedEnergy);

    private:
      bool m_isNormalised;
      float m_displacement;
      float m_DCA;
      float m_trackStubLength;
      float m_nuVertexAvSeparation;
      float m_nuVertexChargeAsymmetry;
      float m_foundConnectionPathway;
      float m_initialGapSize;
      float m_largestGapSize;
      float m_pathwayLength;
      float m_pathwayScatteringAngle2D;
      float m_nShowerHits;
      float m_foundHitRatio;
      float m_scatterAngle;
      float m_openingAngle;
      float m_nuVertexEnergyAsymmetry;
      float m_nuVertexEnergyWeightedMeanRadialDistance;
      float m_showerStartEnergyAsymmetry;
      float m_showerStartMoliereRadius;
      float m_nAmbiguousViews;
      float m_unaccountedEnergy;
  };

    ShowerVarManager(const fhicl::ParameterSet& pset);
    ~ShowerVarManager();

    bool EvaluateShowerVars(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, ShowerVarManager::ShowerVars &showerVars) const;
    void NormaliseShowerVars(ShowerVars &showerVars) const;
    void Reset(ShowerVarManager::ShowerVars &showerVars) const;

  private:
    void FillDisplacement(const art::Event &evt, const art::Ptr<recob::Shower> &shower, ShowerVarManager::ShowerVars &showerVars) const;
    void FillTrackStub(const art::Event &evt, const art::Ptr<recob::Shower> shower, 
        ShowerVarManager::ShowerVars &showerVars) const;
    void FillConnectionPathwayVars(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        ShowerVarManager::ShowerVars &showerVars) const;
    void FillNuVertexAvSeparation(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        ShowerVarManager::ShowerVars &showerVars) const;
    void FillNuVertexChargeAsymmetry(const art::Event &evt, const art::Ptr<recob::Shower> &shower, 
        ShowerVarManager::ShowerVars &showerVars) const;
    float GetViewNuVertexChargeAsymmetry(const art::Event &evt, const TVector3 &nuVertexPosition, const TVector3 &showerStart, 
        const std::vector<art::Ptr<recob::Hit>> &viewHits, const IvysaurusUtils::PandoraView &pandoraView) const;

    std::string m_recoModuleLabel;
    std::string m_showerModuleLabel;
    std::string m_hitModuleLabel;
    float m_displacementLimit;
    float m_DCALimit;
    float m_trackStubLengthLimit;
    float m_nuVertexAvSeparationLimit;
    float m_nuVertexChargeAsymmetryLimit;
    float m_initialGapSizeLimit;
    float m_largestGapSizeLimit;
    float m_pathwayLengthLimit;
    float m_pathwayScatteringAngle2DLimit;
    float m_nShowerHitsLimit;
    float m_foundHitRatioLimit;
    float m_scatterAngleLimit;
    float m_openingAngleLimit;
    float m_nuVertexEnergyAsymmetryLimit;
    float m_nuVertexEnergyWeightedMeanRadialDistanceLimit;
    float m_showerStartEnergyAsymmetryLimit;
    float m_showerStartMoliereRadiusLimit;
    float m_nAmbiguousViewsLimit;
    float m_unaccountedEnergyLimit;
};

/////////////////////////////////////////////////////////////

inline bool ShowerVarManager::ShowerVars::GetIsNormalised() const
{
    return m_isNormalised;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetDisplacement() const
{
    return m_displacement;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetDCA() const
{
    return m_DCA;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetTrackStubLength() const
{
    return m_trackStubLength;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetNuVertexAvSeparation() const
{
    return m_nuVertexAvSeparation;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetNuVertexChargeAsymmetry() const
{
    return m_nuVertexChargeAsymmetry;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetFoundConnectionPathway() const
{
    return m_foundConnectionPathway;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetInitialGapSize() const
{
    return m_initialGapSize;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetLargestGapSize() const
{
    return m_largestGapSize;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetPathwayLength() const
{
    return m_pathwayLength;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetPathwayScatteringAngle2D() const
{
    return m_pathwayScatteringAngle2D;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetNShowerHits() const
{
    return m_nShowerHits;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetFoundHitRatio() const
{
    return m_foundHitRatio;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetScatterAngle() const
{
    return m_scatterAngle;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetOpeningAngle() const
{
    return m_openingAngle;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetNuVertexEnergyAsymmetry() const
{
    return m_nuVertexEnergyAsymmetry;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetNuVertexEnergyWeightedMeanRadialDistance() const
{
    return m_nuVertexEnergyWeightedMeanRadialDistance;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetShowerStartEnergyAsymmetry() const
{
    return m_showerStartEnergyAsymmetry;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetShowerStartMoliereRadius() const
{
    return m_showerStartMoliereRadius;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetNAmbiguousViews() const
{
    return m_nAmbiguousViews;
}

/////////////////////////////////////////////////////////////

inline float ShowerVarManager::ShowerVars::GetUnaccountedEnergy() const
{
    return m_unaccountedEnergy;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetDisplacement(const float displacement) 
{
    m_displacement = displacement;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetDCA(const float dca)
{
    m_DCA = dca;
}


/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetIsNormalised(const bool isNormalised)
{
    m_isNormalised = isNormalised;
}


/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetTrackStubLength(const float trackStubLength)
{
    m_trackStubLength = trackStubLength;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetNuVertexAvSeparation(const float nuVertexAvSeparation)
{
    m_nuVertexAvSeparation = nuVertexAvSeparation;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetNuVertexChargeAsymmetry(const float nuVertexChargeAsymmetry)
{
    m_nuVertexChargeAsymmetry = nuVertexChargeAsymmetry;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetFoundConnectionPathway(const float foundConnectionPathway) 
{
    m_foundConnectionPathway = foundConnectionPathway;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetInitialGapSize(const float initialGapSize) 
{
    m_initialGapSize = initialGapSize;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetLargestGapSize(const float largestGapSize) 
{
    m_largestGapSize = largestGapSize;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetPathwayLength(const float pathwayLength) 
{
    m_pathwayLength = pathwayLength;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetPathwayScatteringAngle2D(const float pathwayScatteringAngle2D) 
{
    m_pathwayScatteringAngle2D = pathwayScatteringAngle2D;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetNShowerHits(const float nShowerHits) 
{
    m_nShowerHits = nShowerHits;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetFoundHitRatio(const float foundHitRatio) 
{
    m_foundHitRatio =  foundHitRatio;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetScatterAngle(const float scatterAngle) 
{
    m_scatterAngle = scatterAngle;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetOpeningAngle(const float openingAngle) 
{
    m_openingAngle = openingAngle;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetNuVertexEnergyAsymmetry(const float nuVertexEnergyAsymmetry) 
{
    m_nuVertexEnergyAsymmetry = nuVertexEnergyAsymmetry;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetNuVertexEnergyWeightedMeanRadialDistance(const float nuVertexEnergyWeightedMeanRadialDistance) 
{
    m_nuVertexEnergyWeightedMeanRadialDistance = nuVertexEnergyWeightedMeanRadialDistance;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetShowerStartEnergyAsymmetry(const float showerStartEnergyAsymmetry) 
{
    m_showerStartEnergyAsymmetry = showerStartEnergyAsymmetry;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetShowerStartMoliereRadius(const float showerStartMoliereRadius) 
{
    m_showerStartMoliereRadius = showerStartMoliereRadius;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetNAmbiguousViews(const float nAmbiguousViews) 
{
    m_nAmbiguousViews = nAmbiguousViews;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetUnaccountedEnergy(const float unaccountedEnergy) 
{
    m_unaccountedEnergy = unaccountedEnergy;
}

/////////////////////////////////////////////////////////////



}

#endif  // SHOWERVARMANAGER_H
