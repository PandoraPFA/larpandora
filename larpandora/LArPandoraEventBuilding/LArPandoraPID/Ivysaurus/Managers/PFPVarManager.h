/**
 *  @file  larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/PFPVarManager.h
 *
 *  @brief A class to manage the Ivysaurus pfp variables input 
 */
#ifndef PFP_VAR_MANAGER_H
#define PFP_VAR_MANAGER_H
// ART
#include "art/Framework/Principal/Event.h"
// LArSoft
#include "lardataobj/RecoBase/PFParticle.h"

namespace ivysaurus
{
/**
 *  @brief  PFPVarManager class
 */
class PFPVarManager
{
  public:
      /**
       *  @brief  PFPVars class
       */
      class PFPVars
      {
      public:
          /**
           *  @brief  Default constructor
           */          
          PFPVars();

          /**
           *  @brief  Get the number of 2D hits (value, is set)
           */
          float GetN2DHits() const;

          /**
           *  @brief  Get the track/shower score (value, is set)
           */          
          float GetTrackShowerScore() const;

          /**
           *  @brief  Return whether the pfp vars have been normalised
           */
          bool GetIsNormalised() const;

          /**
           *  @brief  Set whether the pfp vars have been normalised
           *
           *  @param whether the pfp vars have been normalised
           */
          void SetIsNormalised(const bool isNormalised);

          /**
           *  @brief  Set the number of 2D hits
           *
           *  @param the number of 2D hits           
           */
          void SetN2DHits(const float n2DHits);

            /**
             *  @brief  Set the track/shower score
             *
             *  @param the track/shower score
             */
          void SetTrackShowerScore(const float trackShowerScore);

    private:
          bool m_isNormalised;      ///< whether the pfp vars have been normalised
          float m_n2DHits;          ///< total number of 2D hits
          float m_trackShowerScore; ///< track/shower score (determined by pandora)
  };

    PFPVarManager(const fhicl::ParameterSet& pset);
    ~PFPVarManager();

    /**
     *  @brief  Root function to call to fill the pfp variables for a given PFParticle
     *
     *  @param evt the art event
     *  @param pfparticle the input PFParticle
     *  @param[out] pfpVars the pfp variables to fill
     */    
    void EvaluatePFPVars(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, PFPVarManager::PFPVars &pfpVars) const;

    /**
     *  @brief  Normalise the pfp variables
     *
     *  @param[out] pfpVars the pfp variables to normalise
     */    
    void NormalisePFPVars(PFPVarManager::PFPVars &pfpVars) const;

  private:
    /**
     *  @brief Obtain and fill the track/shower score variable
     *
     *  @param evt the art event
     *  @param pfparticle the input PFParticle
     *  @param pfpVars[out] pfpVars the pfp variables container
     */
    void FillTrackShowerScore(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, PFPVarManager::PFPVars &pfpVars) const;

    /**
     *  @brief Obtain and fill the number of 2D hits variable
     *
     *  @param evt the art event
     *  @param pfparticle the input PFParticle
     *  @param pfpVars[out] pfpVars the pfp variables container
     */    
    void FillN2DHits(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, PFPVarManager::PFPVars &pfpVars) const;

    /**
     *  @brief Normalise an input pfp variable
     *
     *  @param inputPFPVar the input pfp variable
     *  @param mean the mean of the variable (calculated offline)
     *  @param std the standard deviation of the variable (calculated offline)
     *
     *  @return the normalised value
     */    
    float NormalisePFPVar(const float &inputPFPVar, const float mean, const float std) const;    

    std::string m_recoModuleLabel; ///< the pandora label
    float m_n2DHitsMean;           ///< mean for the n2Hits normalisation (calc offline)
    float m_n2DHitsStd;            ///< std for the n2Hits normalisation (calc offline)
    float m_trackShowerScoreMean;  ///< mean for the trackShowerScore normalisation (calc offline)
    float m_trackShowerScoreStd;   ///< std for the trackShowerScore normalisation (calc offline)
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool PFPVarManager::PFPVars::GetIsNormalised() const
{
    return m_isNormalised;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PFPVarManager::PFPVars::GetN2DHits() const
{
    return m_n2DHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PFPVarManager::PFPVars::GetTrackShowerScore() const
{
    return m_trackShowerScore;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void PFPVarManager::PFPVars::SetIsNormalised(const bool isNormalised)
{
    m_isNormalised = isNormalised;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void PFPVarManager::PFPVars::SetN2DHits(const float n2DHits)
{
    m_n2DHits = n2DHits;
}    

//------------------------------------------------------------------------------------------------------------------------------------------

inline void PFPVarManager::PFPVars::SetTrackShowerScore(const float trackShowerScore)
{
    m_trackShowerScore = trackShowerScore;
}    

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace ivysaurus

#endif  // PFP_VAR_MANAGER_H

