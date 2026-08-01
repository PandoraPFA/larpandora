#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/LArPandoraPIDData/PandoraPIDResult.h"

namespace pandorapid
{
    
PandoraPIDResult::PandoraPIDResult()
{
    for (int i = 0; i < static_cast<int>(ParticleType::END); ++i)
    {
        auto particleType = static_cast<ParticleType>(i);
        m_ivysaurusScores.emplace(particleType, -1.f);
    }
}    

//------------------------------------------------------------------------------------------------------------------------------------------

std::pair<ParticleType, float> PandoraPIDResult::GetPredictedParticleType(const ParticleTypeScores &pidScores) const
{
    std::pair<ParticleType, float> prediction = std::make_pair(ParticleType::END, -1.f);
    
    for (const auto &entry : pidScores)
    {
        if ((entry.second > prediction.second) &&
            (std::fabs(entry.second - prediction.second) > std::numeric_limits<float>::epsilon()))              
        {
            prediction = entry;
        }
    }
    
    return prediction;
}

//------------------------------------------------------------------------------------------------------------------------------------------    

void PandoraPIDResult::SetIvysaurusScores(const float muonScore, const float protonScore, const float pionScore,
    const float electronScore, const float photonScore)
{
    m_ivysaurusScores[ParticleType::MUON] = muonScore;
    m_ivysaurusScores[ParticleType::PROTON] = protonScore;
    m_ivysaurusScores[ParticleType::PION] = pionScore;
    m_ivysaurusScores[ParticleType::ELECTRON] = electronScore;
    m_ivysaurusScores[ParticleType::PHOTON] = photonScore;       
}
    
} //namespace pandorapid

