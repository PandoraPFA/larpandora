#ifndef PANDORA_PID_RESULT_H
#define PANDORA_PID_RESULT_H

namespace pandorapid
{
  class PandoraPIDResult
  {
  public:
    double m_ivysaurusMuonScore;
    double m_ivysaurusProtonScore;
    double m_ivysaurusPionScore;
    double m_ivysaurusElectronScore;
    double m_ivysaurusPhotonScore;
    int m_ivysaurusParticleType;
  };
}

#endif
