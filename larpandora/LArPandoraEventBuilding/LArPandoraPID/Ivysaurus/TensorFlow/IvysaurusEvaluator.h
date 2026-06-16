////////////////////////////////////////////////////////////////////////////////////////////////////
//// Class:       IvysaurusEvaluator
//// Authors:     R.Sulej (Robert.Sulej@cern.ch), from DUNE, FNAL/NCBJ, Sept. 2017
///               P.Plonski,                      from DUNE, WUT, Sept. 2017
////              S. Alonso Monsalve,             from DUNE, CERN, Aug. 2018
//// Iterface to run Tensorflow graph saved to a file. First attempts, almost functional.
////
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef IVYSAURUSEVALUATOR_H
#define IVYSAURUSEVALUATOR_H

#include <memory>
#include <vector>
#include <string>

#include "tensorflow/cc/saved_model/loader.h"
#include "tensorflow/cc/saved_model/tag_constants.h"
#include "tensorflow/core/public/session.h"
#include "tensorflow/core/public/session_options.h"
#include "tensorflow/core/framework/logging.h" 

#include "art/Framework/Principal/Event.h"
#include "lardataobj/RecoBase/PFParticle.h"

#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/GridManager.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/TrackVarManager.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/ShowerVarManager.h"

namespace ivysaurus
{

class IvysaurusEvaluator
{
public:
    struct IvysaurusScores
    {
        float m_muonScore;
        float m_protonScore;
        float m_pionScore;
        float m_electronScore;
        float m_photonScore;
        int m_particleType;

        IvysaurusScores() : m_muonScore(-1.f), m_protonScore(-1.f), m_pionScore(-1.f), m_electronScore(-1.f), m_photonScore(-1.f), m_particleType(-1) {};
    };

    IvysaurusEvaluator(fhicl::ParameterSet const &pset);

    IvysaurusScores IvysaurusUseEvaluate(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle);

private:
    std::map<IvysaurusUtils::PandoraView, tensorflow::Tensor> ObtainInputGridTensorMap(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        const bool isStart);

    tensorflow::Tensor ObtainInputTrackTensor(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle);

    tensorflow::Tensor ObtainInputShowerTensor(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle);

    double GetTrackShowerScore(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle);

    std::string m_networkDirectory;
    GridManager m_gridManager;
    TrackVarManager m_trackVarManager;
    ShowerVarManager m_showerVarManager;
    tensorflow::SavedModelBundleLite m_savedModelBundle;
    std::string m_recoModuleLabel;
    int m_nTrackVars;
    int m_nShowerVars;

};

} // namespace ivysaurus

#endif // IVYSAURUSEVALUATOR_H
