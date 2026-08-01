/**
 *  @file  larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/TrackVarManager.h
 *
 *  @brief A class to manage the Ivysaurus track variables input 
 */
#ifndef TRACK_VAR_MANAGER_H
#define TRACK_VAR_MANAGER_H
// ART
#include "art/Framework/Principal/Event.h"
// LArSoft
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

namespace ivysaurus
{
/**
 *  @brief  TrackVarManager class
 */
class TrackVarManager
{
  public:
      /**
       *  @brief  TrackVars class
       */    
      class TrackVars
      {
      public:
          /**
           *  @brief  Default constructor
           */
          TrackVars();

          /**
           *  @brief  Get the number of track children (value, is set)
           */
          std::pair<float, bool> GetNTrackChildren() const;

          /**
           *  @brief  Get the number of shower children (value, is set)
           */
          std::pair<float, bool> GetNShowerChildren() const;

          /**
           *  @brief  Get the number of grand children (value, is set)
           */
          std::pair<float, bool> GetNGrandChildren() const;

          /**
           *  @brief  Get the number of hits of the highest hit child (value, is set)
           */
          std::pair<float, bool> GetNChildHits() const;

          /**
           *  @brief  Get the energy of the highest hit child (value, is set)
           */          
          std::pair<float, bool> GetChildEnergy() const;

          /**
           *  @brief  Get the track score of the highest hit child (value, is set)
           */             
          std::pair<float, bool> GetChildTrackScore() const;

          /**
           *  @brief  Get the track length (value, is set)
           */
          std::pair<float, bool> GetTrackLength() const;

          /**
           *  @brief  Get the std of the defelection along the track (value, is set)
           */
          std::pair<float, bool> GetWobble() const;

          /**
           *  @brief  Get the fabs(byRange - byMCS)/byMCS comparison (value, is set)
           */
          std::pair<float, bool> GetMomentumComparison() const;

          /**
           *  @brief  Return whether the track vars have been normalised
           *
           *  @return whether the track vars have been normalised
           */          
          bool GetIsNormalised() const;

          /**
           *  @brief  Set whether the track vars have been normalised
           *
           *  @param whether the track vars have been normalised
           */
          void SetIsNormalised(const bool isNormalised);

          /**
           *  @brief  Set the number of track children
           *
           *  @param the number of track children
           */
          void SetNTrackChildren(const float nTrackChildren);

          /**
           *  @brief  Set the number of shower children
           *
           *  @param the number of shower children
           */
          void SetNShowerChildren(const float nShowerChildren);

          /**
           *  @brief  Set the number of grand children
           *
           *  @param the number of grand children
           */
          void SetNGrandChildren(const float nGrandChildren);

          /**
           *  @brief  Set the number of hits of the highest hit child
           *
           *  @param the number of hits of the highest hit child
           */
          void SetNChildHits(const float nChildHits);

          /**
           *  @brief  Set the energy of the highest hit child
           *
           *  @param the energy of the highest hit child
           */
          void SetChildEnergy(const float childEnergy);

          /**
           *  @brief  Set the track score of the highest hit child
           *
           *  @param the track score of the highest hit child
           */
          void SetChildTrackScore(const float trackScore);

          /**
           *  @brief  Set the track length
           *
           *  @param the track length
           */
          void SetTrackLength(const float trackLength);

          /**
           *  @brief  Set the std of the defelection along the track (wobble)
           *
           *  @param the std of the defelection along the track (wobble)
           */
          void SetWobble(const float wobble);

          /**
           *  @brief  Set the fabs(byRange - byMCS)/byMCS comparison value
           *
           *  @param the fabs(byRange - byMCS)/byMCS comparison value
           */
          void SetMomentumComparison(const float momComparison);
    private:
          bool m_isNormalised;                         ///< whether the track vars have been normalised
          std::pair<float, bool> m_nTrackChildren;     ///< the number of track children (value, is set)
          std::pair<float, bool> m_nShowerChildren;    ///< the number of shower children (value, is set)
          std::pair<float, bool> m_nGrandChildren;     ///< the number of grand children (value, is set)
          std::pair<float, bool> m_nChildHits;         ///< the number of hits of the highest hit child (value, is set)
          std::pair<float, bool> m_childEnergy;        ///< the energy of the highest hit child (value, is set)
          std::pair<float, bool> m_childTrackScore;    ///< the track score of the highest hit child (value, is set)
          std::pair<float, bool> m_trackLength;        ///< the track length (value, is set)
          std::pair<float, bool> m_wobble;             ///< the std of the defelection along the track (value, is set)
          std::pair<float, bool> m_momentumComparison; ///< the fabs(byRange - byMCS)/byMCS comparison (value, is set)
  };

    TrackVarManager(const fhicl::ParameterSet& pset);
    ~TrackVarManager();

    /**
     *  @brief  Root function to call to fill the track variables for a given PFParticle
     *
     *  @param evt the art event
     *  @param pfparticle the input PFParticle
     *  @param[out] trackVars the track variables to fill
     *
     *  @return whether the track variables could be evaluated (sometimes track fits fail)
     */
    bool EvaluateTrackVars(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, TrackVarManager::TrackVars &trackVars) const;

    /**
     *  @brief  Normalise the track variables
     *
     *  @param[out] trackVars the track variables to normalise
     */
    void NormaliseTrackVars(TrackVarManager::TrackVars &trackVars) const;

  private:
    /**
     *  @brief  Fill the hierarchy related track variables
     *
     *  @param evt the art event
     *  @param pfparticle the input PFParticle
     *  @param[out] trackVars the track variables container
     */
    void FillHierarchyInfo(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, TrackVarManager::TrackVars &trackVars) const;

    /**
     *  @brief  Get the corrected energy of an input PFParticle
     *          by summing hit charge and applying lifetime and recombination corrections
     *
     *  @param evt the art event
     *  @param pfparticle the input PFParticle
     *
     *  @return the corrected energy of the input PFParticle
     */
    float GetPFPEnergy(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle) const;

    /**
     *  @brief Get the track/shower score of an input PFParticle
     *
     *  @param evt the art event
     *  @param pfparticle the input PFParticle
     *
     *  @return the track/shower score of the input PFParticle
     */
    float GetTrackScore(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle) const;

    /**
     *  @brief Fill the track length variable
     *
     *  @param track the input track
     *  @param[out] trackVars the output track variable container
     */
    void FillTrackLength(const art::Ptr<recob::Track> &track, TrackVarManager::TrackVars &trackVars) const;

    /**
     *  @brief Calculate (and fill) the std of the deviation (the opening angle between the directions of successive trajectory points) along the track
     *
     *  @param track the input track
     *  @param[ou] trackVars the track variable container to fill
     */    
    void FillWobble(const art::Ptr<recob::Track> &track, TrackVarManager::TrackVars &trackVars) const;

    /**
     *  @brief Calculate the similarity of the range based and MCS based momentum estimators (assuming muon)
     *         and fill the corresponding trackVars variable
     *
     *  @param track the input track
     *  @param [out] trackVars the track variable container to fill
     */
    void FillTrackMomentum(const art::Ptr<recob::Track> &track, TrackVarManager::TrackVars &trackVars) const;

    /**
     *  @brief Normalise an input track variable, throw an error if the variable has not been filled
     *
     *  @param inputTrackVar the input track variable (variable, has been filled)
     *  @param mean the mean of the variable (calculated offline)
     *  @param std the standard deviation of the variable (calculated offline)
     *
     *  @return the normalised value
     */
    float NormaliseTrackVar(const std::pair<float, bool> &inputTrackVar, const float mean, const float std) const;    

    std::string m_recoModuleLabel;         ///< the pandora label
    std::string m_trackModuleLabel;        ///< the track producer label
    int m_planeIDForEnergy;                ///< ID of plane used in energy-based calculations    
    float m_recombFactor;                  ///< the recombination factor
    calo::CalorimetryAlg m_calorimetryAlg; ///< the calorimetry algorithm used to correct hit energy
    float m_minTrackLengthMCS;             ///< the minimum track length required for the MCS calculation
    float m_maxTrackLengthMCS;             ///< the maximum track length required for the MCS calculation
    float m_intTrkMomRange;                ///< the intercept for the y=mx+c correction of the range-based estimated momentum
    float m_gradTrkMomRange;               ///< the gradient for the y=mx+c correction of the range-based estimated momentum
    float m_intTrkMomMCS;                  ///< the intercept for the y=mx+c correction of the MCS-based estimated momentum
    float m_gradTrkMomMCS;                 ///< the gradient for the y=mx+c correction of the range-based estimated momentum
    float m_nTrackChildrenMean;            ///< mean for the nTrackChildren normalisation (calc offline)
    float m_nTrackChildrenStd;             ///< std for the nTrackChildren normalisation (calc offline)
    float m_nShowerChildrenMean;           ///< mean for the nShowerChildren normalisation (calc offline)
    float m_nShowerChildrenStd;            ///< std for the nShowerChildren normalisation (calc offline)
    float m_nGrandChildrenMean;            ///< mean for the nGrandChildren normalisation (calc offline)
    float m_nGrandChildrenStd;             ///< std for the nGrandChildren normalisation (calc offline)
    float m_nChildHitsMean;                ///< mean for the nChildHits normalisation (calc offline)
    float m_nChildHitsStd;                 ///< std for the nChildHits normalisation (calc offline)
    float m_childEnergyMean;               ///< mean for the childEnergy normalisation (calc offline)
    float m_childEnergyStd;                ///< std for the childenergy normalisation (calc offline)
    float m_childTrackScoreMean;           ///< mean for the childTrackScore normalisation (calc offline)
    float m_childTrackScoreStd;            ///< std for the childTrackScore normalisation (calc offline)
    float m_trackLengthMean;               ///< mean for the trackLength normalisation (calc offline)
    float m_trackLengthStd;                ///< std for the trackLength normalisation (calc offline) 
    float m_wobbleMean;                    ///< mean for the wobble normalisation (calc offline)
    float m_wobbleStd;                     ///< std for the wobble normalisation (calc offline)
    float m_momentumComparisonMean;        ///< mean for the momComparison normalisation (calc offline)
    float m_momentumComparisonStd;         ///< std for the momComparison normalisation (calc offline)
};

//------------------------------------------------------------------------------------------------------------------------------------------    

inline bool TrackVarManager::TrackVars::GetIsNormalised() const
{
    return m_isNormalised;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
inline std::pair<float, bool> TrackVarManager::TrackVars::GetNTrackChildren() const
{
    return m_nTrackChildren;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
inline std::pair<float, bool> TrackVarManager::TrackVars::GetNShowerChildren() const
{
    return m_nShowerChildren;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
inline std::pair<float, bool> TrackVarManager::TrackVars::GetNGrandChildren() const
{
    return m_nGrandChildren;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::pair<float, bool> TrackVarManager::TrackVars::GetNChildHits() const
{
    return m_nChildHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::pair<float, bool> TrackVarManager::TrackVars::GetChildEnergy() const
{
    return m_childEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::pair<float, bool> TrackVarManager::TrackVars::GetChildTrackScore() const
{
    return m_childTrackScore;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::pair<float, bool> TrackVarManager::TrackVars::GetTrackLength() const
{
    return m_trackLength;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::pair<float, bool> TrackVarManager::TrackVars::GetWobble() const
{
    return m_wobble;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::pair<float, bool> TrackVarManager::TrackVars::GetMomentumComparison() const
{
    return m_momentumComparison;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackVarManager::TrackVars::SetIsNormalised(const bool isNormalised)
{
    m_isNormalised = isNormalised;
}

//------------------------------------------------------------------------------------------------------------------------------------------    

inline void TrackVarManager::TrackVars::SetNTrackChildren(const float nTrackChildren)
{
    m_nTrackChildren.first = nTrackChildren;
    m_nTrackChildren.second = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackVarManager::TrackVars::SetNShowerChildren(const float nShowerChildren)
{
    m_nShowerChildren.first = nShowerChildren;
    m_nShowerChildren.second = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackVarManager::TrackVars::SetNGrandChildren(const float nGrandChildren)
{
    m_nGrandChildren.first = nGrandChildren;
    m_nGrandChildren.second = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackVarManager::TrackVars::SetNChildHits(const float nChildHits)
{
    m_nChildHits.first = nChildHits;
    m_nChildHits.second = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackVarManager::TrackVars::SetChildEnergy(const float childEnergy)
{
    m_childEnergy.first = childEnergy;
    m_childEnergy.second = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackVarManager::TrackVars::SetChildTrackScore(const float childTrackScore)
{
    m_childTrackScore.first = childTrackScore;
    m_childTrackScore.second = true;
}    

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackVarManager::TrackVars::SetTrackLength(const float trackLength)
{
    m_trackLength.first = trackLength;
    m_trackLength.second = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackVarManager::TrackVars::SetWobble(const float wobble)
{
    m_wobble.first = wobble;
    m_wobble.second = true;    
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackVarManager::TrackVars::SetMomentumComparison(const float momComparison)
{
    m_momentumComparison.first = momComparison;
    m_momentumComparison.second = true;
}

} // namespace ivysaurus

#endif  // TRACK_VAR_MANAGER_H

