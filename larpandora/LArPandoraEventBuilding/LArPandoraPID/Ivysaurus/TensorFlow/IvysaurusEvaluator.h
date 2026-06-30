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

#include <torch/torch.h>

#include "art/Framework/Principal/Event.h"
#include "lardataobj/RecoBase/PFParticle.h"

#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/GridManager.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/PFPVarManager.h"
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
    //void ObtainGridMask(const GridManager::Grid &grid, std::vector<std::vector<bool>> &mask) const;
    
    std::map<IvysaurusUtils::PandoraView, torch::Tensor> ObtainInputGridTensorMap(const std::map<IvysaurusUtils::PandoraView, GridManager::Grid> &gridMap);
    
    std::map<IvysaurusUtils::PandoraView, torch::Tensor> ObtainGridMaskMap(const std::map<IvysaurusUtils::PandoraView, GridManager::Grid> &gridMap);

    torch::Tensor ObtainInputTrackTensor(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle);

    torch::Tensor ObtainInputShowerTensor(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle);

    bool IsContained(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle);
    
    std::string m_containedNetName;
    std::string m_exitingNetName;    
    GridManager m_gridManager;
    PFPVarManager m_pfpVarManager;
    TrackVarManager m_trackVarManager;
    ShowerVarManager m_showerVarManager;
    std::string m_recoModuleLabel;
    std::string m_trackModuleLabel;
    std::string m_showerModuleLabel;     
    int m_nTrackVars;
    int m_nShowerVars;
    float m_fvMinX;
    float m_fvMaxX;
    float m_fvMinY;
    float m_fvMaxY;
    float m_fvMinZ;
    float m_fvMaxZ;

    torch::jit::script::Module m_containedModel;
    torch::jit::script::Module m_exitingModel;    
};

} // namespace ivysaurus

#endif // IVYSAURUSEVALUATOR_H
