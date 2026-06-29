////////////////////////////////////////////////////////////////////////
/// \file    PFPVarManager.h
/// \brief   A class to manage the Ivysaurus track variables input 
/// \author  Isobel Mawby - i.mawby1@lancaster.ac.uk
////////////////////////////////////////////////////////////////////////

#ifndef PFPMANAGER_H
#define PFPMANAGER_H

#include <vector>
#include <string>
#include <map>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/PFParticle.h"

namespace ivysaurus
{

class PFPVarManager
{
  public:
  class PFPVars
  {
    public:
      PFPVars();

      float GetN2DHits() const;
      float GetTrackShowerScore() const;
      bool GetIsNormalised() const;

      void SetIsNormalised(const bool isNormalised);
      void SetN2DHits(const float n2DHits);
      void SetTrackShowerScore(const float trackShowerScore);

    private:
      bool m_isNormalised;
      float m_n2DHits;
      float m_trackShowerScore;
  };

    PFPVarManager(const fhicl::ParameterSet& pset);
    ~PFPVarManager();

    bool EvaluatePFPVars(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, PFPVarManager::PFPVars &pfpVars) const;
    void NormalisePFPVars(PFPVarManager::PFPVars &pfpVars) const;

  private:
    void FillTrackShowerScore(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, PFPVarManager::PFPVars &pfpVars) const;
    void FillN2DHits(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, PFPVarManager::PFPVars &pfpVars) const;
    float NormalisePFPVar(const float &inputPFPVar, const float mean, const float std) const;    

    std::string m_recoModuleLabel;
    float m_n2DHitsMean;
    float m_n2DHitsStd;    
    float m_trackShowerScoreMean;
    float m_trackShowerScoreStd;
};

/////////////////////////////////////////////////////////////

inline bool PFPVarManager::PFPVars::GetIsNormalised() const
{
    return m_isNormalised;
}

/////////////////////////////////////////////////////////////

inline float PFPVarManager::PFPVars::GetN2DHits() const
{
    return m_n2DHits;
}

/////////////////////////////////////////////////////////////

inline float PFPVarManager::PFPVars::GetTrackShowerScore() const
{
    return m_trackShowerScore;
}

/////////////////////////////////////////////////////////////

inline void PFPVarManager::PFPVars::SetIsNormalised(const bool isNormalised)
{
    m_isNormalised = isNormalised;
}

/////////////////////////////////////////////////////////////

inline void PFPVarManager::PFPVars::SetN2DHits(const float n2DHits)
{
    m_n2DHits = n2DHits;
}    

/////////////////////////////////////////////////////////////

inline void PFPVarManager::PFPVars::SetTrackShowerScore(const float trackShowerScore)
{
    m_trackShowerScore = trackShowerScore;
}    

/////////////////////////////////////////////////////////////

}

#endif  // PFPVARMANAGER_H
