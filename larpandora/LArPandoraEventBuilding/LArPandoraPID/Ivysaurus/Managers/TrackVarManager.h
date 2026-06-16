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

      float GetNTrackChildren() const;
      float GetNShowerChildren() const;
      float GetNGrandChildren() const;
      float GetNChildHits() const;
      float GetChildEnergy() const;
      float GetChildTrackScore() const;
      float GetTrackLength() const;
      float GetWobble() const;
      float GetMomentumComparison() const;
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
      float m_nTrackChildren;
      float m_nShowerChildren;
      float m_nGrandChildren;
      float m_nChildHits;
      float m_childEnergy;
      float m_childTrackScore;
      float m_trackLength;
      float m_wobble;
      float m_momentumComparison;
  };

    TrackVarManager(const fhicl::ParameterSet& pset);
    ~TrackVarManager();

    bool EvaluateTrackVars(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, TrackVarManager::TrackVars &trackVars) const;
    void NormaliseTrackVars(TrackVarManager::TrackVars &trackVars) const;

    void Reset(TrackVarManager::TrackVars &trackVars) const;

  private:
    void FillHierarchyInfo(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, TrackVarManager::TrackVars &trackVars) const;
    float GetChildEnergy(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle) const;
    float GetChildTrackScore(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle) const;
    void FillTrackLength(const art::Ptr<recob::Track> &track, TrackVarManager::TrackVars &trackVars) const;
    void FillWobble(const art::Ptr<recob::Track> &track, TrackVarManager::TrackVars &trackVars) const;
    void FillTrackMomentum(const art::Ptr<recob::Track> &track, TrackVarManager::TrackVars &trackVars) const;

    std::string m_recoModuleLabel;
    std::string m_trackModuleLabel;
    float m_childSeparation;
    float m_recombFactor;
    calo::CalorimetryAlg m_calorimetryAlg;
    float m_minTrackLengthMCS;
    float m_maxTrackLengthMCS;
    float m_intTrkMomRange;
    float m_gradTrkMomRange;
    float m_intTrkMomMCS;
    float m_gradTrkMomMCS;
    float m_nTrackChildrenLimit;
    float m_nShowerChildrenLimit;
    float m_nGrandChildrenLimit;
    float m_nChildHitsLimit;
    float m_childEnergyLimit;
    float m_childTrackScoreLimit;
    float m_trackLengthLimit;
    float m_wobbleLimit;
    float m_momentumComparisonLimit;
};

/////////////////////////////////////////////////////////////

inline bool TrackVarManager::TrackVars::GetIsNormalised() const
{
    return m_isNormalised;
}

/////////////////////////////////////////////////////////////

inline float TrackVarManager::TrackVars::GetNTrackChildren() const
{
    return m_nTrackChildren;
}

/////////////////////////////////////////////////////////////

inline float TrackVarManager::TrackVars::GetNShowerChildren() const
{
    return m_nShowerChildren;
}

/////////////////////////////////////////////////////////////

inline float TrackVarManager::TrackVars::GetNGrandChildren() const
{
    return m_nGrandChildren;
}

/////////////////////////////////////////////////////////////

inline float TrackVarManager::TrackVars::GetNChildHits() const
{
    return m_nChildHits;
}

/////////////////////////////////////////////////////////////

inline float TrackVarManager::TrackVars::GetChildEnergy() const
{
    return m_childEnergy;
}

/////////////////////////////////////////////////////////////

inline float TrackVarManager::TrackVars::GetChildTrackScore() const
{
    return m_childTrackScore;
}

/////////////////////////////////////////////////////////////

inline float TrackVarManager::TrackVars::GetTrackLength() const
{
    return m_trackLength;
}

/////////////////////////////////////////////////////////////

inline float TrackVarManager::TrackVars::GetWobble() const
{
    return m_wobble;
}

/////////////////////////////////////////////////////////////

inline float TrackVarManager::TrackVars::GetMomentumComparison() const
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
    m_nTrackChildren = nTrackChildren;
}

/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetNShowerChildren(const float nShowerChildren)
{
    m_nShowerChildren = nShowerChildren;
}

/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetNGrandChildren(const float nGrandChildren)
{
    m_nGrandChildren = nGrandChildren;
}

/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetNChildHits(const float nChildHits)
{
    m_nChildHits = nChildHits;
}

/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetChildEnergy(const float childEnergy)
{
    m_childEnergy = childEnergy;
}

/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetChildTrackScore(const float trackScore)
{
    m_childTrackScore = trackScore;
}

/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetTrackLength(const float trackLength)
{
    m_trackLength = trackLength;
}

/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetWobble(const float wobble)
{
    m_wobble = wobble;
}

/////////////////////////////////////////////////////////////

inline void TrackVarManager::TrackVars::SetMomentumComparison(const float momComparison)
{
    m_momentumComparison = momComparison;
}




}

#endif  // TRACKVARMANAGER_H
