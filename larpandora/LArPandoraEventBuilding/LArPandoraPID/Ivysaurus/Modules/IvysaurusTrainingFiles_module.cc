/**
 *  @file   dunereco/Ivysaurus/Modules/IvysaurusTrainingFiles_module.cc
 *
 *  @brief  This module uses the analysis utilities to demonstrate 
 *          some of their usage. This can be used as a basis for 
 *          writing analysis code using these tools
 */

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "TTree.h"
#include "TVector3.h"

#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/GridManager.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/PFPVarManager.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/TrackVarManager.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/ShowerVarManager.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Utils/IvysaurusUtils.h"

#include <fstream>
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace ivysaurus
{

/**
 *  @brief  IvysaurusTrainingFiles class
 */
class IvysaurusTrainingFiles : public art::EDAnalyzer
{
public:
   IvysaurusTrainingFiles(fhicl::ParameterSet const &pset);
   virtual ~IvysaurusTrainingFiles();

   void beginJob();
   void endJob();
   void analyze(const art::Event &evt);
   std::vector<int> GetDeltaRays(const art::Event &evt);

   void Reset();

private:

  // Trees
  TTree *m_tree;

  // Tree variables
  int m_run;
  int m_subrun;
  int m_event;
  // General PFP
  int m_truePDG;
  float m_completeness;
  float m_purity;
  int m_nSpacePoints;
  float m_endX;
  float m_endY;
  float m_endZ;
  bool m_isPrimary;
  // Plotting
  std::vector<std::vector<double>> m_spacePoints;
  std::vector<std::vector<double>> m_projectionsU;
  std::vector<std::vector<double>> m_projectionsV;
  std::vector<std::vector<double>> m_projectionsW;
  std::vector<float> m_startDriftBoundariesU;
  std::vector<float> m_startDriftBoundariesV;
  std::vector<float> m_startDriftBoundariesW;
  std::vector<float> m_endDriftBoundariesU;
  std::vector<float> m_endDriftBoundariesV;
  std::vector<float> m_endDriftBoundariesW;
  std::vector<float> m_startWireBoundariesU;
  std::vector<float> m_startWireBoundariesV;
  std::vector<float> m_startWireBoundariesW;
  std::vector<float> m_endWireBoundariesU;
  std::vector<float> m_endWireBoundariesV;
  std::vector<float> m_endWireBoundariesW;
  std::vector<std::vector<float>> m_startGridValuesU;
  std::vector<std::vector<float>> m_startGridValuesV;
  std::vector<std::vector<float>> m_startGridValuesW;
  std::vector<std::vector<float>> m_endGridValuesU;
  std::vector<std::vector<float>> m_endGridValuesV;
  std::vector<std::vector<float>> m_endGridValuesW;
  //PFPVars
  float m_pfpN2DHits;
  float m_pfpTrackShowerScore;
  // TrackVars
  int m_trackVarsSuccessful;
  float m_nTrackChildren;
  float m_nShowerChildren;
  float m_nGrandChildren;
  float m_nChildHits;
  float m_childEnergy;
  float m_childTrackScore;
  float m_trackLength;
  float m_trackWobble;
  float m_trackMomComparison;
  // ShowerVars
  int m_showerVarsSuccessful;
  float m_showerDisplacement;
  float m_DCA;
  float m_trackStubLength;
  float m_nuVertexAvSeparation;
  float m_nuVertexChargeAsymmetry;
  // Managers
  GridManager m_gridManager;
  PFPVarManager m_pfpVarManager;
  TrackVarManager m_trackVarManager;
  ShowerVarManager m_showerVarManager;
  // FCL module labels
  std::string m_hitModuleLabel;
  std::string m_recoModuleLabel;
  std::string m_trackModuleLabel;
  std::string m_showerModuleLabel;    
  // Module variables
  float m_completenessThreshold;
  float m_purityThreshold;
  bool m_writeVisualisationInfo;
};

DEFINE_ART_MODULE(IvysaurusTrainingFiles)

} // namespace ivysaurus

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "larpandora/LArPandoraUtils/PandoraEventUtils.h"
#include "larpandora/LArPandoraUtils/PandoraHitUtils.h"
#include "larpandora/LArPandoraUtils/PandoraPFParticleUtils.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Utils/IvysaurusUtils.h"

#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"

#include <iostream>
#include <random>

namespace ivysaurus
{

IvysaurusTrainingFiles::IvysaurusTrainingFiles(fhicl::ParameterSet const &pset) : 
    art::EDAnalyzer(pset),
    m_gridManager(pset.get<fhicl::ParameterSet>("GridManager")),
    m_pfpVarManager(pset.get<fhicl::ParameterSet>("PFPVarManager")),
    m_trackVarManager(pset.get<fhicl::ParameterSet>("TrackVarManager")),
    m_showerVarManager(pset.get<fhicl::ParameterSet>("ShowerVarManager")),
    m_hitModuleLabel(pset.get<std::string>("HitModuleLabel")),
    m_recoModuleLabel(pset.get<std::string>("RecoModuleLabel")),    
    m_trackModuleLabel(pset.get<std::string>("TrackModuleLabel")),
    m_showerModuleLabel(pset.get<std::string>("ShowerModuleLabel")),      
    m_completenessThreshold(pset.get<float>("CompletenessThreshold")),
    m_purityThreshold(pset.get<float>("PurityThreshold")),
    m_writeVisualisationInfo(pset.get<bool>("WriteVisualisationInfo"))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

IvysaurusTrainingFiles::~IvysaurusTrainingFiles()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void IvysaurusTrainingFiles::analyze(const art::Event &evt)
{
    art::ServiceHandle<cheat::ParticleInventoryService> piServ;
    const std::vector<art::Ptr<recob::PFParticle>> pfparticles = lar_pandora::PandoraEventUtils::GetPFParticles(evt, m_recoModuleLabel);

    // Get the neutrino PFP
    if (!lar_pandora::PandoraEventUtils::HasNeutrino(evt, m_recoModuleLabel))
        return;
    
    art::Ptr<recob::PFParticle> nuPFP = lar_pandora::PandoraEventUtils::GetNeutrino(evt, m_recoModuleLabel);    
    const std::vector<art::Ptr<recob::PFParticle>> &nuChildPFPs = lar_pandora::PandoraPFParticleUtils::GetChildParticles(nuPFP, evt, m_recoModuleLabel);

    for (const art::Ptr<recob::PFParticle> &pfparticle : pfparticles)
    {
        // Make some checks first
        if (!lar_pandora::PandoraPFParticleUtils::HasTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel) &&
            !lar_pandora::PandoraPFParticleUtils::HasShower(pfparticle, evt, m_recoModuleLabel, m_showerModuleLabel))
        {
            continue;
        }
        if (!lar_pandora::PandoraPFParticleUtils::HasVertex(pfparticle, evt, m_recoModuleLabel)) { continue; }
            
        const std::vector<art::Ptr<recob::SpacePoint>> spacepoints = lar_pandora::PandoraPFParticleUtils::GetSpacePoints(pfparticle, evt, m_recoModuleLabel);
        if (spacepoints.empty()) { continue; }
        
        Reset();

        m_run = evt.run();
        m_subrun = evt.subRun();
        m_event = evt.event();
        m_isPrimary = (std::find(nuChildPFPs.begin(), nuChildPFPs.end(), pfparticle) != nuChildPFPs.end());

        ////////////////////////////////////////////
        // First, let's get the truth information...
        ////////////////////////////////////////////  
        const std::vector<art::Ptr<recob::Hit>> pfpHits = lar_pandora::PandoraPFParticleUtils::GetHits(pfparticle, evt, m_recoModuleLabel);
        const std::vector<art::Ptr<recob::Hit>> eventHitList = lar_pandora::PandoraEventUtils::GetHits(evt, m_hitModuleLabel);

        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
        const int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfpHits, 1);

        if (TruthMatchUtils::Valid(g4id))
        {
            // If it isn't a PDG that we care about, move on...
            m_truePDG = piServ->ParticleList().at(g4id)->PdgCode();
            const int absPDG = std::abs(m_truePDG);

            if ((absPDG != 13) && (absPDG != 2212) && (absPDG != 211) && (absPDG != 11) && (absPDG != 22) && (absPDG != 321))
                continue;

            m_completeness = IvysaurusUtils::CompletenessFromTrueParticleID(clockData, pfpHits, eventHitList, g4id);
            m_purity = IvysaurusUtils::HitPurityFromTrueParticleID(clockData, pfpHits, g4id);
        }
        else
        {
            continue;
        }

        ////////////////////////////////////////////
        // Apply truth quality cuts
        ////////////////////////////////////////////  
        if ((m_completeness < m_completenessThreshold) || (m_purity < m_purityThreshold))
            continue;


        ////////////////////////////////////////////
        // Skip if we think it is a DR
        //////////////////////////////////////////// 
        std::vector<int> candidateDeltas = this->GetDeltaRays(evt);
        if (std::find(candidateDeltas.begin(), candidateDeltas.end(), pfparticle->Self()) != candidateDeltas.end())
        {
            continue;
        }
        
        ////////////////////////////////////////////
        // Now, get the endpoints
        ////////////////////////////////////////////
        if (lar_pandora::PandoraPFParticleUtils::HasTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel))
        {
            const art::Ptr<recob::Track> track = lar_pandora::PandoraPFParticleUtils::GetTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel);           
            m_endX = track->End().X(); m_endY = track->End().Y(); m_endZ = track->End().Z();            
        }
      
        ////////////////////////////////////////////
        // Now, get space points into file.. 
        ////////////////////////////////////////////
        m_nSpacePoints = spacepoints.size();

        ////////////////////////////////////////////
        // Initialise grids...
        ////////////////////////////////////////////
        // This is just creating the grid and filling the hit list.
        GridManager::GridMap gridMapStart = m_gridManager.ObtainGridMap(evt, pfparticle, true);

        if (gridMapStart.size() != 3)
            continue;

        GridManager::GridMap gridMapEnd = m_gridManager.ObtainGridMap(evt, pfparticle, false);

        if (gridMapEnd.size() != 3)
            continue;

        ////////////////////////////////////////////
        // Set vars
        ////////////////////////////////////////////
        for (IvysaurusUtils::PandoraView pandoraView : {IvysaurusUtils::PandoraView::TPC_VIEW_U,
             IvysaurusUtils::PandoraView::TPC_VIEW_V, IvysaurusUtils::PandoraView::TPC_VIEW_W})
        {
            const GridManager::Grid &startGrid = gridMapStart.at(pandoraView);
            const GridManager::Grid &endGrid = gridMapEnd.at(pandoraView);
            std::vector<std::vector<float>> &startGridValues = pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_U ? m_startGridValuesU : 
                pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_V ? m_startGridValuesV : m_startGridValuesW;
            std::vector<std::vector<float>> &endGridValues = pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_U ? m_endGridValuesU : 
                pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_V ? m_endGridValuesV : m_endGridValuesW;
            
            const unsigned int dimensions = startGrid.GetAxisDimensions();
            startGridValues = std::vector<std::vector<float>>(dimensions, std::vector<float>(dimensions, 0.f));
            endGridValues = std::vector<std::vector<float>>(dimensions, std::vector<float>(dimensions, 0.f));            

            for (unsigned int driftIndex = 0; driftIndex < dimensions; ++driftIndex)
            {
                for (unsigned int wireIndex = 0; wireIndex <  dimensions; ++wireIndex)
                {
                    startGridValues[driftIndex][wireIndex] = startGrid.GetGridValues().at(driftIndex).at(wireIndex).first;
                    endGridValues[driftIndex][wireIndex] = endGrid.GetGridValues().at(driftIndex).at(wireIndex).first;
                }
            }
        }

        ////////////////////////////////////////////
        // Now fill the pfp variables
        //////////////////////////////////////////// 
        PFPVarManager::PFPVars pfpVars;
        m_pfpVarManager.EvaluatePFPVars(evt, pfparticle, pfpVars);
        m_pfpN2DHits = pfpVars.GetN2DHits();
        m_pfpTrackShowerScore = pfpVars.GetTrackShowerScore();

        ////////////////////////////////////////////
        // Now fill the track variables
        ////////////////////////////////////////////  
        TrackVarManager::TrackVars trackVars;
        m_trackVarsSuccessful = m_trackVarManager.EvaluateTrackVars(evt, pfparticle, trackVars) ? 1 : 0;
        m_nTrackChildren = trackVars.GetNTrackChildren().first;
        m_nShowerChildren = trackVars.GetNShowerChildren().first;
        m_nGrandChildren = trackVars.GetNGrandChildren().first;
        m_nChildHits = trackVars.GetNChildHits().first;
        m_childEnergy = trackVars.GetChildEnergy().first;
        m_childTrackScore = trackVars.GetChildTrackScore().first;
        m_trackLength = trackVars.GetTrackLength().first;
        m_trackWobble = trackVars.GetWobble().first;
        m_trackMomComparison = trackVars.GetMomentumComparison().first;

        ////////////////////////////////////////////
        // Now fill the shower variables
        ////////////////////////////////////////////  
        ShowerVarManager::ShowerVars showerVars;
        m_showerVarsSuccessful = m_showerVarManager.EvaluateShowerVars(evt, pfparticle, showerVars) ? 1 : 0;
        m_showerDisplacement = showerVars.GetDisplacement().first;
        m_DCA = showerVars.GetDCA().first;
        m_trackStubLength = showerVars.GetTrackStubLength().first;
        m_nuVertexAvSeparation = showerVars.GetNuVertexAvSeparation().first;
        m_nuVertexChargeAsymmetry = showerVars.GetNuVertexChargeAsymmetry().first;


        ////////////////////////////////////////////
        // Finally, plotting stuff
        ////////////////////////////////////////////
        if (m_writeVisualisationInfo)
        {
            // Grid boundaries
            for (IvysaurusUtils::PandoraView pandoraView : {IvysaurusUtils::PandoraView::TPC_VIEW_U,
                IvysaurusUtils::PandoraView::TPC_VIEW_V, IvysaurusUtils::PandoraView::TPC_VIEW_W})
            {
                const GridManager::Grid &startGrid = gridMapStart.at(pandoraView);
                const GridManager::Grid &endGrid = gridMapEnd.at(pandoraView);           
                std::vector<float> &startDriftBoundaries = pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_U ? m_startDriftBoundariesU :
                    pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_V ? m_startDriftBoundariesV : m_startDriftBoundariesW;
                std::vector<float> &endDriftBoundaries = pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_U ? m_endDriftBoundariesU :
                    pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_V ? m_endDriftBoundariesV : m_endDriftBoundariesW;
                std::vector<float> &startWireBoundaries = pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_U ? m_startWireBoundariesU : 
                    pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_V ? m_startWireBoundariesV : m_startWireBoundariesW;
                std::vector<float> &endWireBoundaries = pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_U ? m_endWireBoundariesU : 
                    pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_V ? m_endWireBoundariesV : m_endWireBoundariesW;
            
                startDriftBoundaries = startGrid.GetDriftBoundaries();
                endDriftBoundaries = endGrid.GetDriftBoundaries();
                startWireBoundaries = startGrid.GetWireBoundaries();
                endWireBoundaries = endGrid.GetWireBoundaries();
            }
            // Hit info
            std::vector<art::Ptr<recob::Hit>> allHits = lar_pandora::PandoraPFParticleUtils::GetHits(pfparticle, evt, m_recoModuleLabel);        
            const std::vector<art::Ptr<recob::PFParticle>> pfpChildren = lar_pandora::PandoraPFParticleUtils::GetChildParticles(pfparticle, evt, m_recoModuleLabel);
            for (const art::Ptr<recob::PFParticle> &childPFP : pfpChildren)
            {
                const std::vector<art::Ptr<recob::Hit>> &childHits = lar_pandora::PandoraPFParticleUtils::GetHits(childPFP, evt, m_recoModuleLabel);
                allHits.insert(allHits.end(), childHits.begin(), childHits.end());
            }
            for (const art::Ptr<recob::Hit> hit : allHits)
            {
                const IvysaurusUtils::PandoraView thisPandoraView = IvysaurusUtils::GetPandora2DView(hit);
                std::vector<std::vector<double>> &projections = thisPandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_U ? m_projectionsU :
                    thisPandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_V ? m_projectionsV : m_projectionsW;

                // Get 2D hit position
                float width = 0.f;
                TVector3 pandoraHitPosition = TVector3(0.f, 0.f, 0.f);
                IvysaurusUtils::ObtainPandoraHitPositionAndWidth(evt, hit, thisPandoraView, pandoraHitPosition, width);
                projections.push_back({pandoraHitPosition.X(), pandoraHitPosition.Y(), pandoraHitPosition.Z()});
            }
            // Spacepoints
            for (art::Ptr<recob::SpacePoint> spacepoint : spacepoints)
                m_spacePoints.push_back({spacepoint->XYZ()[0], spacepoint->XYZ()[1], spacepoint->XYZ()[2]});
        }             

        m_tree->Fill();
    }
}


//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<int> IvysaurusTrainingFiles::GetDeltaRays(const art::Event &evt)
{
    art::ServiceHandle<cheat::ParticleInventoryService> piServ;
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);

    // Find MCParticle 'muon' counts
    std::vector<int> deltaRays;
    std::map<int, std::vector<art::Ptr<recob::PFParticle>>> muonTrackIDCounts;
    const std::vector<art::Ptr<recob::PFParticle>> pfparticles = lar_pandora::PandoraEventUtils::GetPFParticles(evt, m_recoModuleLabel);
    for (const art::Ptr<recob::PFParticle> &pfparticle : pfparticles)
    {
        const std::vector<art::Ptr<recob::Hit>> pfpHits = lar_pandora::PandoraPFParticleUtils::GetHits(pfparticle, evt, m_recoModuleLabel);
        const int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfpHits, 1);

        if (TruthMatchUtils::Valid(g4id))
        {
            // If it isn't a PDG that we care about, move on...
            int pdg = std::abs(piServ->ParticleList().at(g4id)->PdgCode());

            if (pdg == 13)
                muonTrackIDCounts[g4id].emplace_back(pfparticle);
        }
    }

    if (muonTrackIDCounts.empty())
        return deltaRays;

    // Search for repeated 'muons' and rule out 'obvious' delta rays
    for (const auto &entry : muonTrackIDCounts)
    {
        if (entry.second.size() < 2)
            continue;

        // Make sure we always keep a 'muon'
        int highestHits(0), highestHitPFP(-1);
        for (const art::Ptr<recob::PFParticle> &deltaCandidate : entry.second)
        {
            const std::vector<art::Ptr<recob::Hit>> pfpHits = lar_pandora::PandoraPFParticleUtils::GetHits(deltaCandidate, evt, m_recoModuleLabel);
            if (int(pfpHits.size()) > highestHits)
            {
                highestHits = pfpHits.size();
                highestHitPFP = deltaCandidate->Self();
            }
        }

        for (const art::Ptr<recob::PFParticle> &deltaCandidate : entry.second)
        {
            if (int(deltaCandidate->Self()) == highestHitPFP)
                continue;
            
            const std::vector<art::Ptr<recob::Hit>> pfpHits = lar_pandora::PandoraPFParticleUtils::GetHits(deltaCandidate, evt, m_recoModuleLabel);

            // Is it small and/or shower-like
            if (int(pfpHits.size()) < 90) // 30 in each view
            {
                deltaRays.emplace_back(deltaCandidate->Self());
            }
            else
            {
                const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata = lar_pandora::PandoraPFParticleUtils::GetMetadata(deltaCandidate, evt, m_recoModuleLabel);
                const auto metaMap = metadata->GetPropertiesMap();

                float trackScore = -1.f;
                if (metaMap.find("TrackScore") != metaMap.end())
                    trackScore = metaMap.at("TrackScore");

                if (trackScore < 0.6f)
                    deltaRays.emplace_back(deltaCandidate->Self());
            }
        }
    }

    return deltaRays;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void IvysaurusTrainingFiles::Reset()
{
  const int defaultInt = -999;
  const float defaultFloat = -999.f;

  // General PFP
  m_run = defaultInt;
  m_subrun = defaultInt;
  m_event = defaultInt;
  m_truePDG = defaultInt;
  m_completeness = defaultFloat;
  m_purity = defaultFloat;
  m_nSpacePoints = defaultInt;
  m_endX = defaultFloat;
  m_endY = defaultFloat;
  m_endZ = defaultFloat;
  m_isPrimary = false;
  // Plotting
  m_spacePoints.clear();
  m_projectionsU.clear();
  m_projectionsV.clear();
  m_projectionsW.clear();
  m_startDriftBoundariesU.clear();
  m_startDriftBoundariesV.clear();
  m_startDriftBoundariesW.clear();
  m_endDriftBoundariesU.clear();
  m_endDriftBoundariesV.clear();
  m_endDriftBoundariesW.clear();
  m_startWireBoundariesU.clear();
  m_startWireBoundariesV.clear();
  m_startWireBoundariesW.clear();
  m_endWireBoundariesU.clear();
  m_endWireBoundariesV.clear();
  m_endWireBoundariesW.clear();
  m_startGridValuesU.clear();
  m_startGridValuesV.clear();
  m_startGridValuesW.clear();
  m_endGridValuesU.clear();
  m_endGridValuesV.clear();
  m_endGridValuesW.clear();
  // PFPVars
  m_pfpN2DHits = defaultFloat;  
  m_pfpTrackShowerScore = defaultFloat;
  // TrackVars
  m_trackVarsSuccessful = 0;
  m_nTrackChildren = defaultFloat;
  m_nShowerChildren = defaultFloat;
  m_nGrandChildren = defaultFloat;
  m_nChildHits = defaultFloat;
  m_childEnergy = defaultFloat;
  m_childTrackScore = defaultFloat;
  m_trackLength = defaultFloat;
  m_trackWobble = defaultFloat;
  m_trackMomComparison = defaultFloat;
  // ShowerVars
  m_showerVarsSuccessful = 0;
  m_showerDisplacement = defaultFloat;
  m_DCA = defaultFloat;
  m_trackStubLength = defaultFloat;
  m_nuVertexAvSeparation = defaultFloat;
  m_nuVertexChargeAsymmetry = defaultFloat;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void IvysaurusTrainingFiles::beginJob()
{
    // Implementation of optional member function here.
    art::ServiceHandle<art::TFileService> tfs;
    m_tree = tfs->make<TTree>("ivysaur", "Ivysaur");
    // General PFP
    m_tree->Branch("Run", &m_run);
    m_tree->Branch("Subrun", &m_subrun);
    m_tree->Branch("Event", &m_event);
    m_tree->Branch("TruePDG", &m_truePDG);
    m_tree->Branch("Completeness", &m_completeness);
    m_tree->Branch("Purity", &m_purity);
    m_tree->Branch("NSpacePoints", &m_nSpacePoints);
    m_tree->Branch("EndX", &m_endX);
    m_tree->Branch("EndY", &m_endY);
    m_tree->Branch("EndZ", &m_endZ);    
    m_tree->Branch("IsPrimary", &m_isPrimary, "IsPrimary/O");
    // Plotting
    m_tree->Branch("SpacePoints", &m_spacePoints);
    m_tree->Branch("ProjectionsU", &m_projectionsU);
    m_tree->Branch("ProjectionsV", &m_projectionsV);
    m_tree->Branch("ProjectionsW", &m_projectionsW);
    m_tree->Branch("StartDriftBoundariesU", &m_startDriftBoundariesU);
    m_tree->Branch("StartDriftBoundariesV", &m_startDriftBoundariesV);
    m_tree->Branch("StartDriftBoundariesW", &m_startDriftBoundariesW);
    m_tree->Branch("EndDriftBoundariesU", &m_endDriftBoundariesU);
    m_tree->Branch("EndDriftBoundariesV", &m_endDriftBoundariesV);
    m_tree->Branch("EndDriftBoundariesW", &m_endDriftBoundariesW);
    m_tree->Branch("StartWireBoundariesU", &m_startWireBoundariesU);
    m_tree->Branch("StartWireBoundariesV", &m_startWireBoundariesV);
    m_tree->Branch("StartWireBoundariesW", &m_startWireBoundariesW);
    m_tree->Branch("EndWireBoundariesU", &m_endWireBoundariesU);
    m_tree->Branch("EndWireBoundariesV", &m_endWireBoundariesV);
    m_tree->Branch("EndWireBoundariesW", &m_endWireBoundariesW);
    m_tree->Branch("StartGridU", &m_startGridValuesU);
    m_tree->Branch("StartGridV", &m_startGridValuesV);
    m_tree->Branch("StartGridW", &m_startGridValuesW);
    m_tree->Branch("EndGridU", &m_endGridValuesU);
    m_tree->Branch("EndGridV", &m_endGridValuesV);
    m_tree->Branch("EndGridW", &m_endGridValuesW);
    // PFPVars
    m_tree->Branch("PFPN2DHits", &m_pfpN2DHits);
    m_tree->Branch("PFPTrackShowerScore", &m_pfpTrackShowerScore);
    // TrackVars
    m_tree->Branch("TrackVarsSuccessful", &m_trackVarsSuccessful);
    m_tree->Branch("NTrackChildren", &m_nTrackChildren);
    m_tree->Branch("NShowerChildren", &m_nShowerChildren);
    m_tree->Branch("NGrandChildren", &m_nGrandChildren);
    m_tree->Branch("NChildHits", &m_nChildHits);
    m_tree->Branch("ChildEnergy", &m_childEnergy);
    m_tree->Branch("ChildTrackScore", &m_childTrackScore);
    m_tree->Branch("TrackLength", &m_trackLength);
    m_tree->Branch("TrackWobble", &m_trackWobble);    
    m_tree->Branch("TrackMomComparison", &m_trackMomComparison);
    // ShowerVars
    m_tree->Branch("ShowerVarsSuccessful", &m_showerVarsSuccessful);
    m_tree->Branch("ShowerDisplacement", &m_showerDisplacement);
    m_tree->Branch("ShowerDCA", &m_DCA);
    m_tree->Branch("ShowerTrackStubLength", &m_trackStubLength);
    m_tree->Branch("ShowerNuVertexAvSeparation", &m_nuVertexAvSeparation);
    m_tree->Branch("ShowerNuVertexChargeAsymmetry", &m_nuVertexChargeAsymmetry);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void IvysaurusTrainingFiles::endJob()
{
}


} //namespace ivysaurus

