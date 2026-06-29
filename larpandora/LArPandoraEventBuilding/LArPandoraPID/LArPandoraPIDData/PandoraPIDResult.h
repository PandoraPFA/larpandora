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

    class PandoraPIDResult
    {
    public:
        PandoraPIDResult();
        std::pair<ParticleType, float> GetPredictedParticleType(const std::map<ParticleType, float> &pidScores) const;
        
        std::map<ParticleType, float> m_ivysaurusScores;
    };

    inline PandoraPIDResult::PandoraPIDResult()
    {
        for (int i = 0; i < static_cast<int>(ParticleType::END); ++i)
         {
             auto particleType = static_cast<ParticleType>(i);
             m_ivysaurusScores.emplace(particleType, -1.f);
         }
    }

    inline std::pair<ParticleType, float> PandoraPIDResult::GetPredictedParticleType(const std::map<ParticleType, float> &pidScores) const
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
}

#endif
