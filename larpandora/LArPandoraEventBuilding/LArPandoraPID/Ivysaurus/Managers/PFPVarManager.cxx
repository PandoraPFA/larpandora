////////////////////////////////////////////////////////////////////////
/// \file    TrackVarManager.cxx
/// \brief   A class to manage the Ivysaurus 2D track variable input 
/// \author  Isobel Mawby - i.mawby1@lancaster.ac.uk
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <random>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"

#include "larpandora/LArPandoraUtils/PandoraPFParticleUtils.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/PFPVarManager.h"

namespace ivysaurus
{

PFPVarManager::PFPVars::PFPVars() :
    m_isNormalised(false),
    m_n2DHits(-1.f),
    m_trackShowerScore(-1.f)
{
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PFPVarManager::PFPVarManager(const fhicl::ParameterSet& pset) :
    m_recoModuleLabel(pset.get<std::string>("RecoModuleLabel")),
    m_n2DHitsMean(pset.get<float>("N2DHitsMean")),
    m_n2DHitsStd(pset.get<float>("N2DHitsStd")),    
    m_trackShowerScoreMean(pset.get<float>("TrackShowerScoreMean")),
    m_trackShowerScoreStd(pset.get<float>("TrackShowerScoreStd"))
{
}

/////////////////////////////////////////////////////////////

PFPVarManager::~PFPVarManager()
{
}

/////////////////////////////////////////////////////////////

bool PFPVarManager::EvaluatePFPVars(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    PFPVarManager::PFPVars &pfpVars) const
{
    FillTrackShowerScore(evt, pfparticle, pfpVars);
    FillN2DHits(evt, pfparticle, pfpVars);

    return true;
}

/////////////////////////////////////////////////////////////

void PFPVarManager::FillTrackShowerScore(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, PFPVarManager::PFPVars &pfpVars) const
{
    const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata = lar_pandora::PandoraPFParticleUtils::GetMetadata(pfparticle, evt, m_recoModuleLabel);
    const auto metaMap = metadata->GetPropertiesMap();

    if (metaMap.find("TrackScore") != metaMap.end())
    {
        const double trackScore = metaMap.at("TrackScore");
        pfpVars.SetTrackShowerScore(trackScore);
    }
}

/////////////////////////////////////////////////////////////

void PFPVarManager::FillN2DHits(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, PFPVarManager::PFPVars &pfpVars) const
{
    const std::vector<art::Ptr<recob::Hit>> pfpHits = lar_pandora::PandoraPFParticleUtils::GetHits(pfparticle, evt, m_recoModuleLabel);
    pfpVars.SetN2DHits(pfpHits.size());
}
 
/////////////////////////////////////////////////////////////

void PFPVarManager::NormalisePFPVars(PFPVarManager::PFPVars &pfpVars) const
{
    if (pfpVars.GetIsNormalised())
        throw cet::exception("ivysaur::PFPVarManager::NormalisePFPVars: ") << "pfp vars are already normalised!";  
    
    pfpVars.SetTrackShowerScore(this->NormalisePFPVar(pfpVars.GetTrackShowerScore(), m_trackShowerScoreMean, m_trackShowerScoreStd));
    pfpVars.SetN2DHits(this->NormalisePFPVar(pfpVars.GetN2DHits(), m_n2DHitsMean, m_n2DHitsStd));        
    pfpVars.SetIsNormalised(true);
}

/////////////////////////////////////////////////////////////    
    
float PFPVarManager::NormalisePFPVar(const float &inputPFPVar, const float mean, const float std) const
{
    const float normalised((inputPFPVar - mean) / std);
    return normalised;
}

/////////////////////////////////////////////////////////////

}
