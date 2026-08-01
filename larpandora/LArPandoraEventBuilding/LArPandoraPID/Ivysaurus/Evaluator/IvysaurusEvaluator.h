/**
 *  @file   larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Evaluator/IvysaurusEvaluator.h
 *
 *  @brief  Class to run the Ivysaurus PID
 */
#ifndef IVYSAURUS_EVALUATOR_H
#define IVYSAURUS_EVALUATOR_H
// ART
#include "art/Framework/Principal/Event.h"
// LArSoft
#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/GridManager.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/PFPVarManager.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/TrackVarManager.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/ShowerVarManager.h"
// ROOT
#include <torch/torch.h>
// C++
#include <memory>
#include <vector>
#include <string>

namespace ivysaurus
{
/**
 *  @brief  IvysaurusEvaluator class
 */
class IvysaurusEvaluator
{
public:
    /**
     *  @brief  IvysaurusScores struct
     */        
    struct IvysaurusScores
    {
        float m_muonScore;     ///< Muon PID score
        float m_protonScore;   ///< Proton PID score
        float m_pionScore;     ///< Pion PID score
        float m_electronScore; ///< Electron PID score
        float m_photonScore;   ///< Photon PID score

        /**
         *  @brief  Default constructor
         */        
        IvysaurusScores() : m_muonScore(-1.f), m_protonScore(-1.f), m_pionScore(-1.f), m_electronScore(-1.f), m_photonScore(-1.f) {};
    };

    /**
     *  @brief  Constructor
     *
     *  @param pset parameter set from ivysaurusevaluator.fcl
     */            
    IvysaurusEvaluator(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Run the Ivysaurus PID on an input PFParticle
     *
     *  @param evt the art event
     *  @param pfparticle the input PFParticle
     *
     *  @return the container of ivysaurus scores 
     */    
    IvysaurusScores IvysaurusUseEvaluate(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle);

private:
    /**
     *  @brief  Convert each grid from an input view->grid map into a torch tensor 
     *
     *  @param gridMap the input view->grid map
     *
     *  @return the output view->grid torch tensor map
     */
    std::map<IvysaurusUtils::PandoraView, torch::Tensor> ObtainInputGridTensorMap(const std::map<IvysaurusUtils::PandoraView, GridManager::Grid> &gridMap);

    /**
     *  @brief  Obtain the grid mask (element filled/unfilled=1/0) torch tensor
     *          for each grid in an input view->grid map
     *
     *  @param gridMap the input view->grid map
     *
     *  @return the output view->mask torch tensor map
     */    
    std::map<IvysaurusUtils::PandoraView, torch::Tensor> ObtainGridMaskMap(const std::map<IvysaurusUtils::PandoraView, GridManager::Grid> &gridMap);

    /**
     *  @brief  Obtain the model input torch tensor of track and pfparticle variables
     *
     *  @param evt the art event
     *  @param pfparticle the PFParticle
     *
     *  @return the model input torch tensor of track and pfparticle variables
     */    
    torch::Tensor ObtainInputTrackTensor(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle);

    /**
     *  @brief  Obtain the model input torch tensor of shower variables
     *
     *  @param evt the art event
     *  @param pfparticle the PFParticle
     *
     *  @return the model input torch tensor of shower variables
     */        
    torch::Tensor ObtainInputShowerTensor(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle);

    /**
     *  @brief  Return whether an input PFParticle is contained(true) or exiting(false).
     *          Boundaries are defined in ivysaurusevaluator.fcl
     *
     *  @param evt the art event
     *  @param pfparticle the PFParticle
     *
     *  @return whether the input PFParticle is contained(true) or exiting(false)
     */    
    bool IsContained(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle);
    
    std::string m_containedNetName;      ///< The name of the network used for contained particles
    std::string m_exitingNetName;        ///< The name of the network used for exiting particles
    GridManager m_gridManager;           ///< The GridManager instance
    PFPVarManager m_pfpVarManager;       ///< The PFPVarManager instance
    TrackVarManager m_trackVarManager;   ///< The TrackVarManager instance
    ShowerVarManager m_showerVarManager; ///< The ShowerVarManager instance
    std::string m_recoModuleLabel;       ///< The Pandora label
    std::string m_trackModuleLabel;      ///< The track producer label
    std::string m_showerModuleLabel;     ///< The shower producer label
    int m_nTrackVars;   ///< The number of track+pfp variables that enter the network
    int m_nShowerVars;  ///< The number of shower variables that enter the network
    float m_fvMinX;     ///< The minimum x-coordinate of the 'is contained' region
    float m_fvMaxX;     ///< The maximum x-coordinate of the 'is contained' region
    float m_fvMinY;     ///< The minimum y-coordinate of the 'is contained' region
    float m_fvMaxY;     ///< The maximum y-coordinate of the 'is contained' region
    float m_fvMinZ;     ///< The minimum z-coordinate of the 'is contained' region
    float m_fvMaxZ;     ///< The maximum z-coordinate of the 'is contained' region

    torch::jit::script::Module m_containedModel; ///< The network used for contained particles
    torch::jit::script::Module m_exitingModel;   ///< The network used for exiting particles
};

} // namespace ivysaurus

#endif // IVYSAURUS_EVALUATOR_H

