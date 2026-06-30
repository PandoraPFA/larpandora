////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       IvysaurusEvaluator
// Authors:     R.Sulej (Robert.Sulej@cern.ch), from DUNE, FNAL/NCBJ, Sept. 2017
//              P.Plonski,                      from DUNE, WUT, Sept. 2017
//              S.Alonso-Monsalve,              from DUNE, CERN, Aug. 2018
// Iterface to run Tensorflow graph saved to a file. First attempts, quite functional.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Principal/Event.h"

#include <torch/script.h>
#include <torch/torch.h>

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "larpandora/LArPandoraUtils/PandoraPFParticleUtils.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/GridManager.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/PFPVarManager.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/TrackVarManager.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/ShowerVarManager.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Utils/IvysaurusUtils.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/TensorFlow/IvysaurusEvaluator.h"

/////////////////////////////////////////////////////////////

namespace ivysaurus
{

IvysaurusEvaluator::IvysaurusEvaluator(fhicl::ParameterSet const &pset) :
    m_containedNetName(pset.get<std::string>("ContainedNetName")),
    m_exitingNetName(pset.get<std::string>("ExitingNetName")),     
    m_gridManager(pset.get<fhicl::ParameterSet>("GridManager")),
    m_pfpVarManager(pset.get<fhicl::ParameterSet>("PFPVarManager")),
    m_trackVarManager(pset.get<fhicl::ParameterSet>("TrackVarManager")),
    m_showerVarManager(pset.get<fhicl::ParameterSet>("ShowerVarManager")),
    m_recoModuleLabel(pset.get<std::string>("RecoModuleLabel")),
    m_trackModuleLabel(pset.get<std::string>("TrackModuleLabel")),
    m_showerModuleLabel(pset.get<std::string>("ShowerModuleLabel")),        
    m_nTrackVars(pset.get<int>("NTrackVars")),
    m_nShowerVars(pset.get<int>("NShowerVars")),
    m_fvMinX(pset.get<float>("FVMinX")),
    m_fvMaxX(pset.get<float>("FVMaxX")),
    m_fvMinY(pset.get<float>("FVMinY")),
    m_fvMaxY(pset.get<float>("FVMaxY")),
    m_fvMinZ(pset.get<float>("FVMinZ")),
    m_fvMaxZ(pset.get<float>("FVMaxZ"))
{  
    try
    {
        std::string containedNetPath, exitingNetPath;
        cet::search_path sP("FW_SEARCH_PATH");
        sP.find_file(m_containedNetName, containedNetPath);
        sP.find_file(m_exitingNetName, exitingNetPath);        
        
        m_containedModel = torch::jit::load(containedNetPath);
        m_exitingModel = torch::jit::load(exitingNetPath);        

        // Set the model to evaluation mode.
        // This should have been done during the model export, but we do it here just in case.
        // This ensures that layers like dropout and batch normalization behave correctly during inference.
        m_containedModel.eval();
        m_exitingModel.eval();        
    }
    catch (const std::exception &e)
    {
        std::cout << "Error loading the TorchScript model \'"  << "\':\n" << e.what() << std::endl;
        return;
    }    
}

/////////////////////////////////////////////////////////////

ivysaurus::IvysaurusEvaluator::IvysaurusScores ivysaurus::IvysaurusEvaluator::IvysaurusUseEvaluate(const art::Event &evt, 
    const art::Ptr<recob::PFParticle> &pfparticle)
{
    IvysaurusScores ivysaurusScores;

    // Make some checks first
    if (!lar_pandora::PandoraPFParticleUtils::HasTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel) &&
        !lar_pandora::PandoraPFParticleUtils::HasShower(pfparticle, evt, m_recoModuleLabel, m_showerModuleLabel))
    {
        return ivysaurusScores;
    }
    const std::vector<art::Ptr<recob::SpacePoint>> spacepoints = lar_pandora::PandoraPFParticleUtils::GetSpacePoints(pfparticle, evt, m_recoModuleLabel);
    if (spacepoints.empty()) { return ivysaurusScores; }

    // Get grid maps
    GridManager::GridMap startGridMap = m_gridManager.ObtainGridMap(evt, pfparticle, true);
    if (startGridMap.size() != 3) { return ivysaurusScores; }
    GridManager::GridMap endGridMap = m_gridManager.ObtainGridMap(evt, pfparticle, false);
    if (endGridMap.size() != 3) { return ivysaurusScores; }    
    
    // Normalise grids
    for (IvysaurusUtils::PandoraView pandoraView : {IvysaurusUtils::PandoraView::TPC_VIEW_U, 
         IvysaurusUtils::PandoraView::TPC_VIEW_V, IvysaurusUtils::PandoraView::TPC_VIEW_W})
    {
        GridManager::Grid &startGrid = startGridMap.at(pandoraView);
        GridManager::Grid &endGrid = endGridMap.at(pandoraView);        
        m_gridManager.NormaliseGrid(startGrid);
        m_gridManager.NormaliseGrid(endGrid);        
    }
    
    // Convert to torch tensors...
    std::map<IvysaurusUtils::PandoraView, torch::Tensor> startGridTensorMap = ObtainInputGridTensorMap(startGridMap);
    std::map<IvysaurusUtils::PandoraView, torch::Tensor> endGridTensorMap = ObtainInputGridTensorMap(endGridMap);
    std::map<IvysaurusUtils::PandoraView, torch::Tensor> startMaskMap = this->ObtainGridMaskMap(startGridMap);
    std::map<IvysaurusUtils::PandoraView, torch::Tensor> endMaskMap = this->ObtainGridMaskMap(endGridMap); 
    
    // Obtain the input track variable tensors (PFPVars included in track)
    torch::Tensor trackVarTensor = ObtainInputTrackTensor(evt, pfparticle);

    // Obtain the input shower variable tensors
    torch::Tensor showerVarTensor = ObtainInputShowerTensor(evt, pfparticle);

    // Run the model
    const bool isContained = this->IsContained(evt, pfparticle);
    torch::jit::script::Module &ivysaurus = isContained ? m_containedModel : m_exitingModel;
    torch::NoGradGuard guard;
    torch::Tensor output = ivysaurus.forward({startGridTensorMap.at(IvysaurusUtils::TPC_VIEW_U), startMaskMap.at(IvysaurusUtils::TPC_VIEW_U),
            endGridTensorMap.at(IvysaurusUtils::TPC_VIEW_U), endMaskMap.at(IvysaurusUtils::TPC_VIEW_U),
            startGridTensorMap.at(IvysaurusUtils::TPC_VIEW_V), startMaskMap.at(IvysaurusUtils::TPC_VIEW_V),
            endGridTensorMap.at(IvysaurusUtils::TPC_VIEW_V), endMaskMap.at(IvysaurusUtils::TPC_VIEW_V),
            startGridTensorMap.at(IvysaurusUtils::TPC_VIEW_W), startMaskMap.at(IvysaurusUtils::TPC_VIEW_W),
            endGridTensorMap.at(IvysaurusUtils::TPC_VIEW_W), endMaskMap.at(IvysaurusUtils::TPC_VIEW_W),
            trackVarTensor, showerVarTensor}).toTensor();
        
    // Remember to softmax!
    torch::Tensor probs = torch::softmax(output, 1);
    ivysaurusScores.m_muonScore = probs[0][0].item<float>();
    ivysaurusScores.m_protonScore = probs[0][1].item<float>();
    ivysaurusScores.m_pionScore = probs[0][2].item<float>();
    ivysaurusScores.m_electronScore = probs[0][3].item<float>();
    ivysaurusScores.m_photonScore = probs[0][4].item<float>();

    float highestScore = -std::numeric_limits<float>::max();
    int count = 0;

    for (float score : {ivysaurusScores.m_muonScore, ivysaurusScores.m_protonScore, ivysaurusScores.m_pionScore, 
        ivysaurusScores.m_electronScore, ivysaurusScores.m_photonScore})
    {
        if (score > highestScore)
        {
            highestScore = score;
            ivysaurusScores.m_particleType = count;
        }

        ++count;
    }

    return ivysaurusScores;
}

/////////////////////////////////////////////////////////////

bool IvysaurusEvaluator::IsContained(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle)
{
    if (!lar_pandora::PandoraPFParticleUtils::HasTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel))
        return false;

    const art::Ptr<recob::Track> track = lar_pandora::PandoraPFParticleUtils::GetTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel);
    const float endX(track->End().X()), endY(track->End().Y()), endZ(track->End().Z());

    if ((endX < m_fvMinX) || (endX > m_fvMaxX))
        return false;

    if ((endY < m_fvMinY) || (endY > m_fvMaxY))
        return false;

    if ((endZ < m_fvMinZ) || (endZ > m_fvMaxZ))
        return false;

    return true;
}

/////////////////////////////////////////////////////////////

std::map<IvysaurusUtils::PandoraView, torch::Tensor> IvysaurusEvaluator::ObtainInputGridTensorMap(const std::map<IvysaurusUtils::PandoraView, GridManager::Grid> &gridMap)
{
    std::map<IvysaurusUtils::PandoraView, torch::Tensor> tensorViewMap;

    for (IvysaurusUtils::PandoraView pandoraView : {IvysaurusUtils::PandoraView::TPC_VIEW_U,
         IvysaurusUtils::PandoraView::TPC_VIEW_V, IvysaurusUtils::PandoraView::TPC_VIEW_W})
    {
        const GridManager::Grid &grid(gridMap.at(pandoraView));
        torch::Tensor gridTensor = torch::zeros({1, grid.GetAxisDimensions(), grid.GetAxisDimensions(), 1});

        for (unsigned int driftIndex = 0; driftIndex < grid.GetAxisDimensions(); ++driftIndex)
        {
            for (unsigned int wireIndex = 0; wireIndex < grid.GetAxisDimensions(); ++wireIndex)
            {
                gridTensor[0][driftIndex][wireIndex][0] = grid.GetGridValues().at(driftIndex).at(wireIndex).first;
            }
        }

        tensorViewMap[pandoraView] = gridTensor;
    }

    return tensorViewMap;
}

/////////////////////////////////////////////////////////////

std::map<IvysaurusUtils::PandoraView, torch::Tensor> IvysaurusEvaluator::ObtainGridMaskMap(const std::map<IvysaurusUtils::PandoraView, GridManager::Grid> &gridMap)
{
    std::map<IvysaurusUtils::PandoraView, torch::Tensor> maskMap;

    for (IvysaurusUtils::PandoraView pandoraView : {IvysaurusUtils::PandoraView::TPC_VIEW_U,
         IvysaurusUtils::PandoraView::TPC_VIEW_V, IvysaurusUtils::PandoraView::TPC_VIEW_W})
    {
        const GridManager::Grid &grid(gridMap.at(pandoraView));
        torch::Tensor maskTensor = torch::zeros({1, grid.GetAxisDimensions(), grid.GetAxisDimensions(), 1});

        for (unsigned int driftIndex = 0; driftIndex < grid.GetAxisDimensions(); ++driftIndex)
        {
            for (unsigned int wireIndex = 0; wireIndex < grid.GetAxisDimensions(); ++wireIndex)
            {
                maskTensor[0][driftIndex][wireIndex][0] = float(grid.GetGridValues().at(driftIndex).at(wireIndex).second);
            }
        }
        
        maskMap[pandoraView] = maskTensor;
    }

    return maskMap;
}

/////////////////////////////////////////////////////////////

torch::Tensor IvysaurusEvaluator::ObtainInputTrackTensor(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle)
{
    torch::Tensor trackVarTensor = torch::zeros({1, m_nTrackVars});

    PFPVarManager::PFPVars pfpVars;
    m_pfpVarManager.EvaluatePFPVars(evt, pfparticle, pfpVars);
    m_pfpVarManager.NormalisePFPVars(pfpVars);
    
    TrackVarManager::TrackVars trackVars;
    m_trackVarManager.EvaluateTrackVars(evt, pfparticle, trackVars);
    m_trackVarManager.NormaliseTrackVars(trackVars);

    // ATTN: Order is important!
    trackVarTensor[0][0] = trackVars.GetNTrackChildren().first;
    trackVarTensor[0][1] = (trackVars.GetNTrackChildren().second ? 1.f : 0.f);
    trackVarTensor[0][2] = trackVars.GetNShowerChildren().first;
    trackVarTensor[0][3] = (trackVars.GetNShowerChildren().second ? 1.f : 0.f);
    trackVarTensor[0][4] = trackVars.GetNGrandChildren().first;
    trackVarTensor[0][5] = (trackVars.GetNGrandChildren().second ? 1.f : 0.f);
    trackVarTensor[0][6] = trackVars.GetNChildHits().first;
    trackVarTensor[0][7] = (trackVars.GetNChildHits().second ? 1.f : 0.f);
    trackVarTensor[0][8] = trackVars.GetChildEnergy().first;
    trackVarTensor[0][9] = (trackVars.GetChildEnergy().second ? 1.f : 0.f); 
    trackVarTensor[0][10] = trackVars.GetChildTrackScore().first;
    trackVarTensor[0][11] = (trackVars.GetChildTrackScore().second ? 1.f : 0.f);  
    trackVarTensor[0][12] = trackVars.GetTrackLength().first;
    trackVarTensor[0][13] = (trackVars.GetTrackLength().second ? 1.f : 0.f);  
    trackVarTensor[0][14] = trackVars.GetWobble().first;
    trackVarTensor[0][15] = (trackVars.GetWobble().second ? 1.f : 0.f); 
    trackVarTensor[0][16] = trackVars.GetMomentumComparison().first;
    trackVarTensor[0][17] = (trackVars.GetMomentumComparison().second ? 1.f : 0.f);
    trackVarTensor[0][18] = pfpVars.GetN2DHits();
    trackVarTensor[0][19] = pfpVars.GetTrackShowerScore();    
    return trackVarTensor;
}

// /////////////////////////////////////////////////////////////

torch::Tensor IvysaurusEvaluator::ObtainInputShowerTensor(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle)
{
    torch::Tensor showerVarTensor = torch::zeros({1, m_nShowerVars});

    ShowerVarManager::ShowerVars showerVars;
    m_showerVarManager.EvaluateShowerVars(evt, pfparticle, showerVars);
    m_showerVarManager.NormaliseShowerVars(showerVars);

    // ATTN: Order is important!
    showerVarTensor[0][0] = showerVars.GetDisplacement().first;
    showerVarTensor[0][1] = (showerVars.GetDisplacement().second ? 1.f : 0.f);
    showerVarTensor[0][2] = showerVars.GetDCA().first;
    showerVarTensor[0][3] = (showerVars.GetDCA().second ? 1.f : 0.f);
    showerVarTensor[0][4] = showerVars.GetTrackStubLength().first;
    showerVarTensor[0][5] = (showerVars.GetTrackStubLength().second ? 1.f : 0.f);
    showerVarTensor[0][6] = showerVars.GetNuVertexAvSeparation().first;
    showerVarTensor[0][7] = (showerVars.GetNuVertexAvSeparation().second ? 1.f : 0.f);
    showerVarTensor[0][8] = showerVars.GetNuVertexChargeAsymmetry().first;
    showerVarTensor[0][9] = (showerVars.GetNuVertexChargeAsymmetry().second ? 1.f : 0.f);

    return showerVarTensor;
}

}
