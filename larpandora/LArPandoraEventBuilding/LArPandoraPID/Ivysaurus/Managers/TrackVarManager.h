////////////////////////////////////////////////////////////////////////
/// \file    TrackVarManager.h
/// \brief   A class to manage the Ivysaurus track variables input 
/// \author  Isobel Mawby - i.mawby1@lancaster.ac.uk
////////////////////////////////////////////////////////////////////////

#ifndef TRACKMANAGER_H
#define TRACKMANAGER_H

#include <vector>
#include <string>
#include <map>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"

namespace ivysaurus
{

class TrackVarManager
{
  public:
  class TrackVars
  {
    public:
      TrackVars();

      std::pair<float, bool> GetNTrackChildren() const;
      std::pair<float, bool> GetNShowerChildren() const;
      std::pair<float, bool> GetNGrandChildren() const;
      std::pair<float, bool> GetNChildHits() const;
      std::pair<float, bool> GetChildEnergy() const;
      std::pair<float, bool> GetChildTrackScore() const;
      std::pair<float, bool> GetTrackLength() const;
      std::pair<float, bool> GetWobble() const;
      std::pair<float, bool> GetMomentumComparison() const;
      bool GetIsNormalised() const;

      void SetIsNormalised(const bool isNormalised);
      void SetNTrackChildren(const float nTrackChildren);
      void SetNShowerChildren(const float nShowerChildren);
      void SetNGrandChildren(const float nGrandChildren);
      void SetNChildHits(const float nChildHits);
      void SetChildEnergy(const float childEnergy);
      void SetTrackLength(const float trackLength);
      void SetChildTrackScore(const float trackScore);
      void SetWobble(const float wobble);
      void SetMomentumComparison(const float momComparison);

    private:
      bool m_isNormalised;
      std::pair<float, bool> m_nTrackChildren;
      std::pair<float, bool> m_nShowerChildren;
      std::pair<float, bool> m_nGrandChildren;
      std::pair<float, bool> m_nChildHits;
      std::pair<float, bool> m_childEnergy;
      std::pair<float, bool> m_childTrackScore;
      std::pair<float, bool> m_trackLength;
      std::pair<float, bool> m_wobble;
      std::pair<float, bool> m_momentumComparison;
  };

    TrackVarManager(const fhicl::ParameterSet& pset);
    ~TrackVarManager();

    bool EvaluateTrackVars(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, TrackVarManager::TrackVars &trackVars) const;
    void NormaliseTrackVars(TrackVarManager::TrackVars &trackVars) const;

  private:
    void FillHierarchyInfo(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, TrackVarManager::TrackVars &trackVars) const;
    float GetChildEnergy(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle) const;
    float GetTrackScore(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle) const;
    void FillTrackLength(const art::Ptr<recob::Track> &track, TrackVarManager::TrackVars &trackVars) const;
    void FillWobble(const art::Ptr<recob::Track> &track, TrackVarManager::TrackVars &trackVars) const;
    void FillTrackMomentum(const art::Ptr<recob::Track> &track, TrackVarManager::TrackVars &trackVars) const;
    float NormaliseTrackVar(const std::pair<float, bool> &inputTrackVar, const float mean, const float std) const;    

    std::string m_recoModuleLabel;
    std::string m_trackModuleLabel;
    float m_recombFactor;
    calo::CalorimetryAlg m_calorimetryAlg;
    float m_minTrackLengthMCS;
    float m_maxTrackLengthMCS;
    float m_intTrkMomRange;
    float m_gradTrkMomRange;
    float m_intTrkMomMCS;
    float m_gradTrkMomMCS;
    float m_nTrackChildrenMean;
    float m_nTrackChildrenStd;    
    float m_nShowerChildrenMean;
    float m_nShowerChildrenStd;    
    float m_nGrandChildrenMean;
    float m_nGrandChildrenStd;    
    float m_nChildHitsMean;
    float m_nChildHitsStd;    
    float m_childEnergyMean;
    float m_childEnergyStd;    
    float m_childTrackScoreMean;
    float m_childTrackScoreStd;    
    float m_trackLengthMean;
    float m_trackLengthStd;    
    float m_wobbleMean;
    float m_wobbleStd;
    float m_momentumComparisonMean;
    float m_momentumComparisonStd;    
};

/////////////////////////////////////////////////////////////

inline bool TrackVarManager::TrackVars::GetIsNormalised() const
{
    return m_isNormalised;
}

/////////////////////////////////////////////////////////////

inline std::pair<float, bool> TrackVarManager::TrackVars::GetNTrackChildren() const
{
    return m_nTrackChildren;
}

/////////////////////////////////////////////////////////////

inline std::pair<float, bool> TrackVarManager::TrackVars::GetNShowerChildren() const
{
    return m_nShowerChildren;
}

/////////////////////////////////////////////////////////////

inline std::pair<float, bool> TrackVarManager::TrackVars::GetNGrandChildren() const
{
    return m_nGrandChildren;
}

/////////////////////////////////////////////////////////////

inline std::pair<float, bool> TrackVarManager::TrackVars::GetNChildHits() const
{
    return m_nChildHits;
}

/////////////////////////////////////////////////////////////

inline std::pair<float, bool> TrackVarManager::TrackVars::GetChildEnergy() const
{
    return m_childEnergy;
}

/////////////////////////////////////////////////////////////

inline std::pair<float, bool> TrackVarManager::TrackVars::GetChildTrackScore() const
{
    return m_childTrackScore;
}

/////////////////////////////////////////////////////////////

inline std::pair<float, bool> TrackVarManager::TrackVars::GetTrackLength() const
{
    return m_trackLength;
}

/////////////////////////////////////////////////////////////

inline std::pair<float, bool> TrackVarManager::TrackVars::GetWobble() const
{
    return m_wobble;
}

/////////////////////////////////////////////////////////////

inline std::pair<float, bool> TrackVarManager::TrackVars::GetMomentumComparison() const
{
    return m_momentumComparison;
}

/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetIsNormalised(const bool isNormalised)
{
    m_isNormalised = isNormalised;
}
    
/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetNTrackChildren(const float nTrackChildren)
{
    m_nTrackChildren.first = nTrackChildren;
    m_nTrackChildren.second = true;
}

/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetNShowerChildren(const float nShowerChildren)
{
    m_nShowerChildren.first = nShowerChildren;
    m_nShowerChildren.second = true;
}

/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetNGrandChildren(const float nGrandChildren)
{
    m_nGrandChildren.first = nGrandChildren;
    m_nGrandChildren.second = true;
}

/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetNChildHits(const float nChildHits)
{
    m_nChildHits.first = nChildHits;
    m_nChildHits.second = true;
}

/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetChildEnergy(const float childEnergy)
{
    m_childEnergy.first = childEnergy;
    m_childEnergy.second = true;
}

/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetChildTrackScore(const float childTrackScore)
{
    m_childTrackScore.first = childTrackScore;
    m_childTrackScore.second = true;
}    

/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetTrackLength(const float trackLength)
{
    m_trackLength.first = trackLength;
    m_trackLength.second = true;
}

/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetWobble(const float wobble)
{
    m_wobble.first = wobble;
    m_wobble.second = true;    
}

/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetMomentumComparison(const float momComparison)
{
    m_momentumComparison.first = momComparison;
    m_momentumComparison.second = true;
}

/////////////////////////////////////////////////////////////

}

#endif  // TRACKVARMANAGER_H
