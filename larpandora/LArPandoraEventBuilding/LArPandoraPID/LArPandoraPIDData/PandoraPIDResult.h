/**
 *  @file  larpandora/LArPandoraEventBuilding/LArPandoraPID/LArPandoraPIDData/PandoraPIDResult.h
 *
 *  @brief A class to hold the scores of various Pandora PID estimators
 */
#ifndef PANDORA_PID_RESULT_H
#define PANDORA_PID_RESULT_H

#include <map>
#include <cmath>
#include <limits>

namespace pandorapid
{
enum ParticleType {
    MUON,
    PROTON,
    PION,
    KAON,
    ELECTRON,
    PHOTON,
    END
};

typedef std::map<ParticleType, float> ParticleTypeScores;
    
/**
 *  @brief  PandoraPIDResult class
 */        
class PandoraPIDResult
{
public:
    /**
     *  @brief  Default constructor 
     */    
    PandoraPIDResult();

    /**
     *  @brief  Given an input ParticleType->score map, return the ParticleType with the highest score
     *
     *  @param  pidScores the input ParticleType->score map
     *
     *  @return the predicted ParticleType alongside its score
     */    
    std::pair<ParticleType, float> GetPredictedParticleType(const ParticleTypeScores &pidScores) const;

    /**
     *  @brief  Get the Ivysaurus PID ParticleType->score map
     *
     *  @return the Ivysaurus PID ParticleType->score map
     */     
    const ParticleTypeScores& GetIvysaurusScores() const;
    
    /**
     *  @brief  Set Ivysaurus PID scores
     *
     *  @param  muonScore the muon score
     *  @param  protonScore the proton score
     *  @param  pionScore the pion score
     *  @param  electronScore the electron score
     *  @param  photonScore the photon score     
     */     
    void SetIvysaurusScores(const float muonScore, const float protonScore, const float pionScore,
        const float electronScore, const float photonScore);

private:    
    ParticleTypeScores m_ivysaurusScores; ///< The ParticleType->score map for the Ivysaurus PID
};

//------------------------------------------------------------------------------------------------------------------------------------------
    
inline const ParticleTypeScores& PandoraPIDResult::GetIvysaurusScores() const
{
    return m_ivysaurusScores;
}

} // namespace pandorapid

#endif

