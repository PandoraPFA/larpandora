/**
 *  @file  larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/ShowerVarManager.h
 *
 *  @brief A class to manage the Ivysaurus shower variables input 
 *
 */
#ifndef SHOWER_VAR_MANAGER_H
#define SHOWER_VAR_MANAGER_H
// ART
#include "art/Framework/Principal/Event.h"
//LArSoft
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Utils/IvysaurusUtils.h"

namespace ivysaurus
{
/**
 *  @brief  ShowerVarManager class
 */
class ShowerVarManager
{
  public:
      /**
       *  @brief  TrackVars class
       */        
      class ShowerVars
      {
      public:
          /**
           *  @brief  Default constructor
           */          
          ShowerVars();
          
          /**
           *  @brief  Return whether the shower vars have been normalised
           *
           *  @return whether the shower vars have been normalised
           */          
          bool GetIsNormalised() const;

          /**
           *  @brief Get the displacement (value, is set)
           */
          std::pair<float, bool> GetDisplacement() const;

          /**
           *  @brief Get the distance of closest approach (value, is set)
           */          
          std::pair<float, bool> GetDCA() const;

          /**
           *  @brief Get the length of the track stub (value, is set)
           */          
          std::pair<float, bool> GetTrackStubLength() const;

          /**
           *  @brief Get the average separation of the PFParticle's spacepoints from a central axis
           *         defined by the particles startpoint and the parent's endpoint (value, is set)
           */                    
          std::pair<float, bool> GetFromParentAvSep() const;

          /**
           *  @brief Get the asymmetry of charge distributed around a central axis
           *         defined by the particles startpoint and the parent's endpoint (value, is set)
           *         calculated in 2D, maximum value across the 2D views is taken
           */
          std::pair<float, bool> GetFromParentChargeAsym() const;

          /**
           *  @brief  Set whether the shower vars have been normalised
           *
           *  @param whether the shower vars have been normalised
           */
          void SetIsNormalised(const bool isNormalised);

          /**
           *  @brief Set the displacement
           *
           *  @param the displacement
           */
          void SetDisplacement(const float displacement);

          /**
           *  @brief Set the distance of closest approach
           *
           *  @param the distance of closest approach
           */
          void SetDCA(const float dca);

          /**
           *  @brief Set the track stub length
           *
           *  @param the track stub length
           */
          void SetTrackStubLength(const float trackStubLength);

          /**
           *  @brief Set the 'from parent' average separation
           *
           *  @param the 'from parent' average separation
           */          
          void SetFromParentAvSep(const float fromParentAvSep);

          /**
           *  @brief Set the 'from parent' charge asymmetry
           *
           *  @param the 'from parent' charge asymmetry
           */
          void SetFromParentChargeAsym(const float fromParentChargeAsym);
      
    private:
          bool m_isNormalised;                           ///< whether the shower variables have been normalised
          std::pair<float, bool> m_displacement;         ///< the displacement from the parent endpoint (value, is set)
          std::pair<float, bool> m_DCA;                  ///< the distance of closest approach to the parent endpoint (value, is set)
          std::pair<float, bool> m_trackStubLength;      ///< the track stub length (value, is set)
          std::pair<float, bool> m_fromParentAvSep;      ///< the parent endpoint aveage separation (value, is set)
          std::pair<float, bool> m_fromParentChargeAsym; ///< the parent charge asymmetry (value, is set)
  };

    ShowerVarManager(const fhicl::ParameterSet& pset);
    ~ShowerVarManager();

    /**
     *  @brief  Root function to call to fill the shower variables for a given PFParticle
     *
     *  @param evt the art event
     *  @param pfparticle the input PFParticle
     *  @param[out] showerVars the shower variables to fill
     *
     *  @return whether the shower variables could be evaluated (sometimes shower fits fail)
     */    
    bool EvaluateShowerVars(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, ShowerVarManager::ShowerVars &showerVars) const;

    /**
     *  @brief  Normalise the shower variables
     *
     *  @param[out] showerVars the shower variables to normalise
     */    
    void NormaliseShowerVars(ShowerVars &showerVars) const;

  private:

    /**
     *  @brief  Identify the endpoint of the parent particle (neutrino vertex is pfparticle is primary)
     *
     *  @param evt the art event
     *  @param pfparticle the input PFParticle
     *  @param[out] parentEnd the endpoint of the parent PFParticle
     *
     *  @return whether the parent endpoint could be identified (sometimes track fits fail)
     */
    bool GetParentEndpoint(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, TVector3 &parentEnd) const;

    /**
     *  @brief Calculate and fill the displacement varaible
     *
     *  @param evt the art event
     *  @param parentEndpoint the endpoint of the parent PFParticle
     *  @param shower the input shower
     *  @param[out] the shower variable container
     */   
    void FillDisplacement(const art::Event &evt, const TVector3 &parentEndpoint, const art::Ptr<recob::Shower> &shower, ShowerVarManager::ShowerVars &showerVars) const;

    /**
     *  @brief Calculate and fill the track stub length varaible
     *
     *  @param evt the art event
     *  @param shower the input shower
     *  @param[out] the shower variable container
     */   
    void FillTrackStub(const art::Event &evt, const art::Ptr<recob::Shower> shower, 
        ShowerVarManager::ShowerVars &showerVars) const;

    /**
     *  @brief Calculate (and fill) the average separation of the PFParticle's spacepoints from a central axis
     *         defined by the particles startpoint and the parent's endpoint
     *
     *  @param evt the art event
     *  @param parentEndpoint the endpoint of the parent PFParticle
     *  @param pfparticle the input PFParticle
     *  @param shower the input shower
     *  @param[out] the shower variable container
     */
     void FillAvSeparation(const art::Event &evt, const TVector3 &parentEndpoint, const art::Ptr<recob::PFParticle> &pfparticle, 
        const art::Ptr<recob::Shower> &shower, ShowerVarManager::ShowerVars &showerVars) const;

    /**
     *  @brief Calculate (and fill) the asymmetry of the charge around a central axis
     *         defined by the particle's startpoint and the parent's endpoint
     *         in each 2D view, maximum value across the 2D views is taken
     *
     *  @param evt the art event
     *  @param parentEndpoint the endpoint of the parent PFParticle
     *  @param shower the input shower
     *  @param[out] the shower variable container
     */
    void FillChargeAsymmetry(const art::Event &evt, const TVector3 &parentEndpoint, const art::Ptr<recob::Shower> &shower, 
        ShowerVarManager::ShowerVars &showerVars) const;

    /**
     *  @brief Calculate the asymmetry of the charge around a central axis
     *         defined by the particle's startpoint and the parent's endpoint
     *         in a 2D view
     *
     *  @param evt the art event
     *  @param parentEndpoint the endpoint of the parent PFParticle
     *  @param showerStart the start of the input PFParticle
     *  @param viewHits the PFParticle hits in the 2D view
     *  @param pandoraView the input 2D view
     *
     *  @return the charge asymmetry in the specified view
     */
    float GetViewChargeAsymmetry(const art::Event &evt, const TVector3 &parentEndpoint, const TVector3 &showerStart, 
        const std::vector<art::Ptr<recob::Hit>> &viewHits, const IvysaurusUtils::PandoraView &pandoraView) const;

    /**
     *  @brief Normalise an input shower variable, throw an error if the variable has not been filled
     *
     *  @param inputShowerVar the input shower variable (variable, has been filled)
     *  @param mean the mean of the variable (calculated offline)
     *  @param std the standard deviation of the variable (calculated offline)
     *
     *  @return the normalised value
     */    
    float NormaliseShowerVar(const std::pair<float, bool> &inputShowerVar, const float mean, const float std) const;    

    std::string m_recoModuleLabel;    ///< the pandora label
    std::string m_trackModuleLabel;   ///< the track producer label
    std::string m_showerModuleLabel;  ///< the shower producer label
    std::string m_hitModuleLabel;     ///< the hit producer label
    int m_planeIDForEnergy;           ///< ID of plane used in energy-based calculations
    float m_displacementMean;         ///< mean for the displacement normalisation (calc offline)
    float m_displacementStd;          ///< std for the displacement normalisation (calc offline)
    float m_DCAMean;                  ///< mean for the DCA normalisation (calc offline)
    float m_DCAStd;                   ///< std for the DCA normalisation (calc offline)
    float m_trackStubLengthMean;      ///< mean for the track stub normalisation (calc offline)
    float m_trackStubLengthStd;       ///< std for the track stub normalisation (calc offline)
    float m_fromParentAvSepMean;      ///< mean for the fromParentAvSep normalisation (calc offline)
    float m_fromParentAvSepStd;       ///< std for the fromParentAvSep normalisation (calc offline)
    float m_fromParentChargeAsymMean; ///< mean for the fromParentChargeAsym normalisation (calc offline)
    float m_fromParentChargeAsymStd;  ///< std for the fromParentChargeAsym normalisation (calc offline)
};

//------------------------------------------------------------------------------------------------------------------------------------------
    
inline bool ShowerVarManager::ShowerVars::GetIsNormalised() const
{
    return m_isNormalised;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::pair<float, bool> ShowerVarManager::ShowerVars::GetDisplacement() const
{
    return m_displacement;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::pair<float, bool> ShowerVarManager::ShowerVars::GetDCA() const
{
    return m_DCA;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::pair<float, bool> ShowerVarManager::ShowerVars::GetTrackStubLength() const
{
    return m_trackStubLength;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::pair<float, bool> ShowerVarManager::ShowerVars::GetFromParentAvSep() const
{
    return m_fromParentAvSep;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::pair<float, bool> ShowerVarManager::ShowerVars::GetFromParentChargeAsym() const
{
    return m_fromParentChargeAsym;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ShowerVarManager::ShowerVars::SetDisplacement(const float displacement) 
{
    m_displacement.first = displacement;
    m_displacement.second = true;    
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ShowerVarManager::ShowerVars::SetDCA(const float dca)
{
    m_DCA.first = dca;
    m_DCA.second = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ShowerVarManager::ShowerVars::SetTrackStubLength(const float trackStubLength)
{
    m_trackStubLength.first = trackStubLength;
    m_trackStubLength.second = true;    
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ShowerVarManager::ShowerVars::SetFromParentAvSep(const float fromParentAvSep)
{
    m_fromParentAvSep.first = fromParentAvSep;
    m_fromParentAvSep.second = true;    
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ShowerVarManager::ShowerVars::SetFromParentChargeAsym(const float fromParentChargeAsym)
{
    m_fromParentChargeAsym.first = fromParentChargeAsym;
    m_fromParentChargeAsym.second = true;    
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

inline void ShowerVarManager::ShowerVars::SetIsNormalised(const bool isNormalised)
{
    m_isNormalised = isNormalised;
}

} // namespace ivysaurus

#endif  // SHOWER_VAR_MANAGER_H

