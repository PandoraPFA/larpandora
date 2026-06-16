////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       IvysaurusEvaluator
// Authors:     R.Sulej (Robert.Sulej@cern.ch), from DUNE, FNAL/NCBJ, Sept. 2017
//              P.Plonski,                      from DUNE, WUT, Sept. 2017
//              S.Alonso-Monsalve,              from DUNE, CERN, Aug. 2018
// Iterface to run Tensorflow graph saved to a file. First attempts, quite functional.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Principal/Event.h"

#include "tensorflow/cc/saved_model/loader.h"
#include "tensorflow/cc/saved_model/tag_constants.h"
#include "tensorflow/core/public/session.h"
#include "tensorflow/core/public/session_options.h"
#include "tensorflow/core/framework/logging.h" 

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "larpandora/LArPandoraUtils/PandoraPFParticleUtils.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/GridManager.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/TrackVarManager.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/ShowerVarManager.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Utils/IvysaurusUtils.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/TensorFlow/IvysaurusEvaluator.h"

/////////////////////////////////////////////////////////////

ivysaurus::IvysaurusEvaluator::IvysaurusEvaluator(fhicl::ParameterSet const &pset) :
    m_networkDirectory(pset.get<std::string>("NetworkDirectory")), 
    m_gridManager(pset.get<fhicl::ParameterSet>("GridManager")),
    m_trackVarManager(pset.get<fhicl::ParameterSet>("TrackVarManager")),
    m_showerVarManager(pset.get<fhicl::ParameterSet>("ShowerVarManager")),
    m_recoModuleLabel(pset.get<std::string>("RecoModuleLabel")),
    m_nTrackVars(pset.get<int>("NTrackVars")),
    m_nShowerVars(pset.get<int>("NShowerVars"))
{
    // Create dummy options.
    tensorflow::SessionOptions sessionOptions;
    tensorflow::RunOptions runOptions;

    std::string directoryPath;
    cet::search_path sP("FW_SEARCH_PATH");
    sP.find_file(m_networkDirectory, directoryPath);

    //std::cout << "Directory path: " << directoryPath << std::endl;

    // Load the model bundle. (this returns a status code)
    const auto loadResult = tensorflow::LoadSavedModel(sessionOptions, runOptions, directoryPath, 
        {tensorflow::kSavedModelTagServe}, &m_savedModelBundle);

    // Check if loading was okay.
    TF_CHECK_OK(loadResult);
}

/////////////////////////////////////////////////////////////

ivysaurus::IvysaurusEvaluator::IvysaurusScores ivysaurus::IvysaurusEvaluator::IvysaurusUseEvaluate(const art::Event &evt, 
    const art::Ptr<recob::PFParticle> &pfparticle)
{
    IvysaurusScores ivysaurusScores;
    
    // Obtain the input grid tensors
    //std::cout << "Making the grid tensors..." << std::endl;
    std::map<IvysaurusUtils::PandoraView, tensorflow::Tensor> startGridTensorMap = ObtainInputGridTensorMap(evt, pfparticle, true);
    if (startGridTensorMap.size() != 3) { return ivysaurusScores; }
    std::map<IvysaurusUtils::PandoraView, tensorflow::Tensor> endGridTensorMap = ObtainInputGridTensorMap(evt, pfparticle, false);
    if (endGridTensorMap.size() != 3) { return ivysaurusScores; }

    // Obtain the input track variable tensors
    //std::cout << "Making the track variable tensors..." << std::endl;
    tensorflow::Tensor trackVarTensor = ObtainInputTrackTensor(evt, pfparticle);

    // auto trackVarTensorMap = trackVarTensor.tensor<float, 2>();
    // for (int i = 0; i < m_nTrackVars; ++i)
    //     std::cout << "trackVarTensorMap(0, i): " << trackVarTensorMap(0, i) << std::endl;

    // Obtain the input shower variable tensors
    //std::cout << "Making the shower variable tensors..." << std::endl;
    tensorflow::Tensor showerVarTensor = ObtainInputShowerTensor(evt, pfparticle);

    // auto showerVarTensorMap = showerVarTensor.tensor<float, 2>();
    // for (int i = 0; i < m_nShowerVars; ++i)
    //     std::cout << "showerVarTensorMap(0, i): " << showerVarTensorMap(0, i) << std::endl;

    // Link the data with some tags so tensorflow know where to put those data entries.
    //std::cout << "Hooking up input/output layers..." << std::endl;
    std::vector<std::pair<std::string, tensorflow::Tensor>> feedInputs = {{"serve_input_1:0", startGridTensorMap.at(IvysaurusUtils::PandoraView::TPC_VIEW_U)}, 
                                                                          {"serve_input_2:0", endGridTensorMap.at(IvysaurusUtils::PandoraView::TPC_VIEW_U)},
                                                                          {"serve_input_3:0", startGridTensorMap.at(IvysaurusUtils::PandoraView::TPC_VIEW_V)}, 
                                                                          {"serve_input_4:0", endGridTensorMap.at(IvysaurusUtils::PandoraView::TPC_VIEW_V)},
                                                                          {"serve_input_5:0", startGridTensorMap.at(IvysaurusUtils::PandoraView::TPC_VIEW_W)},
                                                                          {"serve_input_6:0", endGridTensorMap.at(IvysaurusUtils::PandoraView::TPC_VIEW_W)},
                                                                          {"serve_input_7:0", trackVarTensor},
                                                                          {"serve_input_8:0", showerVarTensor}};
    std::vector<std::string> fetches = { "StatefulPartitionedCall:0" };

    // We need to store the results somewhere.
    std::vector<tensorflow::Tensor> outputs;

    // Let's run the model...
    //std::cout << "Running the model..." << std::endl;
    auto status = m_savedModelBundle.GetSession()->Run(feedInputs, fetches, {}, &outputs);

    // ... and print out it's predictions.
    //std::cout << "What does the model predict?..." << std::endl;
    if (!status.ok()) { return ivysaurusScores; }

    auto output_map = outputs[0].tensor<float, 2>();
    ivysaurusScores.m_muonScore = output_map(0, 0);
    ivysaurusScores.m_protonScore = output_map(0, 1);
    ivysaurusScores.m_pionScore = output_map(0, 2);
    ivysaurusScores.m_electronScore = output_map(0, 3);
    ivysaurusScores.m_photonScore = output_map(0, 4);

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

    std::cout << "--------------------------------" << std::endl;
    std::cout << "found an ivysaurus score!!!" << std::endl;
    std::cout << "muonScore: " << ivysaurusScores.m_muonScore << std::endl;
    std::cout << "protonScore: " << ivysaurusScores.m_protonScore << std::endl;
    std::cout << "pionScore: " << ivysaurusScores.m_pionScore << std::endl;
    std::cout << "electronScore: " << ivysaurusScores.m_electronScore << std::endl;
    std::cout << "photonScore: " << ivysaurusScores.m_photonScore << std::endl;
    std::cout << "particleType: " << ivysaurusScores.m_particleType << std::endl;

    TF_CHECK_OK(status);

    return ivysaurusScores;
}

/////////////////////////////////////////////////////////////

std::map<IvysaurusUtils::PandoraView, tensorflow::Tensor> ivysaurus::IvysaurusEvaluator::ObtainInputGridTensorMap(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    const bool isStart)
{
    // The tensor view map to fill and output
    std::map<IvysaurusUtils::PandoraView, tensorflow::Tensor> tensorViewMap;

    // Fill the grids
    GridManager::GridMap gridMap = m_gridManager.ObtainGridMap(evt, pfparticle, isStart);

    if (gridMap.size() != 3)
        return tensorViewMap;

    m_gridManager.FillGrids(evt, pfparticle, gridMap);

    // Convert to tensors
    for (IvysaurusUtils::PandoraView pandoraView : {IvysaurusUtils::PandoraView::TPC_VIEW_U,
         IvysaurusUtils::PandoraView::TPC_VIEW_V, IvysaurusUtils::PandoraView::TPC_VIEW_W})
    {
        GridManager::Grid &grid(gridMap.at(pandoraView));

        // Create the tensors
        tensorflow::Tensor gridTensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, grid.GetAxisDimensions(), grid.GetAxisDimensions(), 1 }));
        auto tensorMap = gridTensor.tensor<float, 4>();

        for (unsigned int driftIndex = 0; driftIndex < grid.GetAxisDimensions(); ++driftIndex)
        {
            for (unsigned int wireIndex = 0; wireIndex < grid.GetAxisDimensions(); ++wireIndex)
            {
                tensorMap(0, driftIndex, wireIndex, 0) = grid.GetGridValues().at(driftIndex).at(wireIndex);
                //std::cout << "(" << driftIndex << ", " << wireIndex << "): " << grid.GetGridValues().at(driftIndex).at(wireIndex) << std::endl;
            }
        }

        tensorViewMap[pandoraView] = gridTensor;
    }

    return tensorViewMap;
}

/////////////////////////////////////////////////////////////

tensorflow::Tensor ivysaurus::IvysaurusEvaluator::ObtainInputTrackTensor(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle)
{
    tensorflow::Tensor trackVarTensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({1, m_nTrackVars}));
    auto trackVarsTensorMap = trackVarTensor.tensor<float, 2>();
    
    TrackVarManager::TrackVars trackVars;
    m_trackVarManager.EvaluateTrackVars(evt, pfparticle, trackVars);
    m_trackVarManager.NormaliseTrackVars(trackVars);

    // ATTN: Order is important!
    trackVarsTensorMap(0, 0) = trackVars.GetNTrackChildren();
    trackVarsTensorMap(0, 1) = trackVars.GetNShowerChildren();
    trackVarsTensorMap(0, 2) = trackVars.GetNGrandChildren();
    trackVarsTensorMap(0, 3) = trackVars.GetNChildHits();
    trackVarsTensorMap(0, 4) = trackVars.GetChildEnergy();
    trackVarsTensorMap(0, 5) = trackVars.GetChildTrackScore();
    trackVarsTensorMap(0, 6) = trackVars.GetTrackLength();
    trackVarsTensorMap(0, 7) = trackVars.GetWobble();
    trackVarsTensorMap(0, 8) = GetTrackShowerScore(evt, pfparticle);
    trackVarsTensorMap(0, 9) = trackVars.GetMomentumComparison();
    return trackVarTensor;
}

/////////////////////////////////////////////////////////////

tensorflow::Tensor ivysaurus::IvysaurusEvaluator::ObtainInputShowerTensor(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle)
{
    tensorflow::Tensor showerVarTensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({1, m_nShowerVars}));
    auto showerVarsTensorMap = showerVarTensor.tensor<float, 2>();

    ShowerVarManager::ShowerVars showerVars;

    m_showerVarManager.EvaluateShowerVars(evt, pfparticle, showerVars);
    m_showerVarManager.NormaliseShowerVars(showerVars);

    // ATTN: Order is important!
    showerVarsTensorMap(0, 0) = showerVars.GetDisplacement();
    showerVarsTensorMap(0, 1) = showerVars.GetDCA();
    showerVarsTensorMap(0, 2) = showerVars.GetInitialGapSize();
    showerVarsTensorMap(0, 3) = showerVars.GetLargestGapSize();
    showerVarsTensorMap(0, 4) = showerVars.GetPathwayLength();
    showerVarsTensorMap(0, 5) = showerVars.GetPathwayScatteringAngle2D();
    showerVarsTensorMap(0, 6) = showerVars.GetNShowerHits();
    showerVarsTensorMap(0, 7) = showerVars.GetFoundHitRatio();
    showerVarsTensorMap(0, 8) = showerVars.GetScatterAngle();
    showerVarsTensorMap(0, 9) = showerVars.GetOpeningAngle();
    showerVarsTensorMap(0, 10) = showerVars.GetNuVertexEnergyAsymmetry();
    showerVarsTensorMap(0, 11) = showerVars.GetNuVertexEnergyWeightedMeanRadialDistance();
    showerVarsTensorMap(0, 12) = showerVars.GetShowerStartEnergyAsymmetry();
    showerVarsTensorMap(0, 13) = showerVars.GetShowerStartMoliereRadius();
    showerVarsTensorMap(0, 14) = showerVars.GetNAmbiguousViews();
    showerVarsTensorMap(0, 15) = showerVars.GetUnaccountedEnergy();

    return showerVarTensor;
}

/////////////////////////////////////////////////////////////

double ivysaurus::IvysaurusEvaluator::GetTrackShowerScore(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle)
{
    const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata = lar_pandora::PandoraPFParticleUtils::GetMetadata(pfparticle, evt, m_recoModuleLabel);
    const auto metaMap = metadata->GetPropertiesMap();

    if (metaMap.find("TrackScore") == metaMap.end())
        return -1.0;

    const double trackScore = metaMap.at("TrackScore");

    return trackScore;
}

/////////////////////////////////////////////////////////////
