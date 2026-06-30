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
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"

#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

#include "larpandora/LArPandoraUtils/PandoraPFParticleUtils.h"
#include "larpandora/LArPandoraUtils/PandoraHitUtils.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Managers/TrackVarManager.h"

namespace ivysaurus
{

TrackVarManager::TrackVars::TrackVars() : 
        m_isNormalised(false),
        m_nTrackChildren(std::make_pair(-1.f, false)),
        m_nShowerChildren(std::make_pair(-1.f, false)),
        m_nGrandChildren(std::make_pair(-1.f, false)),
        m_nChildHits(std::make_pair(-1.f, false)),
        m_childEnergy(std::make_pair(-1.f, false)),
        m_childTrackScore(std::make_pair(-1.f, false)),
        m_trackLength(std::make_pair(-1.f, false)),
        m_wobble(std::make_pair(-1.f, false)),
        m_momentumComparison(std::make_pair(-1.f, false))
{
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TrackVarManager::TrackVarManager(const fhicl::ParameterSet& pset) :
    m_recoModuleLabel(pset.get<std::string>("RecoModuleLabel")),
    m_trackModuleLabel(pset.get<std::string>("TrackModuleLabel")),
    m_recombFactor(pset.get<float>("RecombFactor")),
    m_calorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
    m_minTrackLengthMCS(pset.get<float>("MinTrackLengthMCS")),
    m_maxTrackLengthMCS(pset.get<float>("MaxTrackLengthMCS")),
    m_intTrkMomRange(pset.get<float>("IntTrkMomRange")),
    m_gradTrkMomRange(pset.get<float>("GradTrkMomRange")),
    m_intTrkMomMCS(pset.get<float>("IntTrkMomMCS")),
    m_gradTrkMomMCS(pset.get<float>("GradTrkMomMCS")),
    m_nTrackChildrenMean(pset.get<float>("NTrackChildrenMean")),
    m_nTrackChildrenStd(pset.get<float>("NTrackChildrenStd")),    
    m_nShowerChildrenMean(pset.get<float>("NShowerChildrenMean")),
    m_nShowerChildrenStd(pset.get<float>("NShowerChildrenStd")),    
    m_nGrandChildrenMean(pset.get<float>("NGrandChildrenMean")),
    m_nGrandChildrenStd(pset.get<float>("NGrandChildrenStd")),    
    m_nChildHitsMean(pset.get<float>("NChildHitsMean")),
    m_nChildHitsStd(pset.get<float>("NChildHitsStd")),    
    m_childEnergyMean(pset.get<float>("ChildEnergyMean")),
    m_childEnergyStd(pset.get<float>("ChildEnergyStd")),    
    m_childTrackScoreMean(pset.get<float>("ChildTrackScoreMean")),
    m_childTrackScoreStd(pset.get<float>("ChildTrackScoreStd")),    
    m_trackLengthMean(pset.get<float>("TrackLengthMean")),
    m_trackLengthStd(pset.get<float>("TrackLengthStd")),    
    m_wobbleMean(pset.get<float>("WobbleMean")),
    m_wobbleStd(pset.get<float>("WobbleStd")),
    m_momentumComparisonMean(pset.get<float>("MomentumComparisonMean")),
    m_momentumComparisonStd(pset.get<float>("MomentumComparisonStd"))
{
}

/////////////////////////////////////////////////////////////

TrackVarManager::~TrackVarManager()
{
}

/////////////////////////////////////////////////////////////

bool TrackVarManager::EvaluateTrackVars(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    TrackVarManager::TrackVars &trackVars) const
{
    if (!lar_pandora::PandoraPFParticleUtils::HasTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel))
        return false;

    const art::Ptr<recob::Track> track = lar_pandora::PandoraPFParticleUtils::GetTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel);

    FillHierarchyInfo(evt, pfparticle, trackVars);
    FillTrackLength(track, trackVars);
    FillWobble(track, trackVars);
    FillTrackMomentum(track, trackVars);

    return true;
}

/////////////////////////////////////////////////////////////

void TrackVarManager::FillHierarchyInfo(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    TrackVarManager::TrackVars &trackVars) const
{
    if (!lar_pandora::PandoraPFParticleUtils::IsTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel))
        return;
    
    const art::Ptr<recob::Track> track = lar_pandora::PandoraPFParticleUtils::GetTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel);
    const TVector3 trackEndpoint = TVector3(track->End().X(), track->End().Y(), track->End().Z());
    const std::vector<art::Ptr<recob::PFParticle>> childPFPs = lar_pandora::PandoraPFParticleUtils::GetChildParticles(pfparticle, evt, m_recoModuleLabel);
    
    int nTracks = 0, nShowers = 0, nGrandChildren = 0;
    int highestHits = -1;
    float highestHitEnergy = -1, highestHitTrackScore = -1; 

    for (const art::Ptr<recob::PFParticle> &childPFP : childPFPs)
    {
        const std::vector<art::Ptr<recob::SpacePoint>> spacepoints = lar_pandora::PandoraPFParticleUtils::GetSpacePoints(childPFP, evt, m_recoModuleLabel);

        if (spacepoints.empty())
            continue;

        if (childPFP->PdgCode() == 13)
            ++nTracks;
        else if (childPFP->PdgCode() == 11)
            ++nShowers;
        
        const std::vector<art::Ptr<recob::PFParticle>> grandChildPFPs = lar_pandora::PandoraPFParticleUtils::GetChildParticles(childPFP, evt, m_recoModuleLabel);
        nGrandChildren += grandChildPFPs.size();

        const std::vector<art::Ptr<recob::Hit>> hits = lar_pandora::PandoraPFParticleUtils::GetHits(childPFP, evt, m_recoModuleLabel);
        if (static_cast<int>(hits.size()) > highestHits)
        {
            highestHits = hits.size();
            highestHitEnergy = GetChildEnergy(evt, pfparticle);
            highestHitTrackScore = this->GetTrackScore(evt, pfparticle);
        }
    }
    
    trackVars.SetNTrackChildren(nTracks);
    trackVars.SetNShowerChildren(nShowers);
    trackVars.SetNGrandChildren(nGrandChildren);

    if (highestHits != -1)
    {
        trackVars.SetNChildHits(highestHits);
        trackVars.SetChildEnergy(highestHitEnergy);
        trackVars.SetChildTrackScore(highestHitTrackScore);
    }
}

/////////////////////////////////////////////////////////////

float TrackVarManager::GetChildEnergy(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle) const
{
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
    const std::vector<art::Ptr<recob::Hit>> &pfpHits = lar_pandora::PandoraPFParticleUtils::GetViewHits(pfparticle, evt, m_recoModuleLabel, 2);
    const double charge(lar_pandora::PandoraHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, pfpHits));
    const double energy = m_calorimetryAlg.ElectronsFromADCArea(charge, 2) / m_recombFactor / util::kGeVToElectrons;

    return energy;
}

/////////////////////////////////////////////////////////////

float TrackVarManager::GetTrackScore(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle) const
{
    const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata = lar_pandora::PandoraPFParticleUtils::GetMetadata(pfparticle, evt, m_recoModuleLabel);
    const auto metaMap = metadata->GetPropertiesMap();

    if (metaMap.find("TrackScore") == metaMap.end())
        return -1.f;

    const float trackScore = metaMap.at("TrackScore");

    return trackScore;
}

/////////////////////////////////////////////////////////////

void TrackVarManager::FillTrackLength(const art::Ptr<recob::Track> &track, TrackVarManager::TrackVars &trackVars) const
{
    const float trackLength = track->Length();
    trackVars.SetTrackLength(trackLength);
}

/////////////////////////////////////////////////////////////
// This is copied from Dom's Pandizzle module
// Follows the MicroBooNE method for MCS    
void TrackVarManager::FillWobble(const art::Ptr<recob::Track> &track, TrackVarManager::TrackVars &trackVars) const
{
  //Get the number of points
  size_t NPoints = track->NumberTrajectoryPoints();
  //Store the directions between adjacent points on a vector
  std::vector<TVector3> directions;
  for (size_t i_point = 0; i_point < NPoints-1; i_point++){
    TVector3 position_i(track->TrajectoryPoint(i_point).position.X(), track->TrajectoryPoint(i_point).position.Y(), track->TrajectoryPoint(i_point).position.Z());
    TVector3 position_iplus1(track->TrajectoryPoint(i_point+1).position.X(), track->TrajectoryPoint(i_point+1).position.Y(), track->TrajectoryPoint(i_point+1).position.Z());
    TVector3 direction = (position_iplus1-position_i).Unit();
    directions.push_back(direction);
  }

  //Loop through the direction and compare adjacent elements
  std::vector<double> deflection_angles;
  for (size_t i_dir = 0; i_dir < directions.size()-1; i_dir++){
    //Aim: rotate both direction so that the first direction is parallel to the z-axis.  
    //Then take the x-projection of scattered track and calculate the angle between that and the first direction
    TVector3 z_axis(0,0,1);
    TVector3 direction_first = directions[i_dir];
    TVector3 direction_second = directions[i_dir+1];

    //Ignore if either direction is 0 (i.e. not 1)
    if (direction_first.Mag() < 0.999 || direction_second.Mag() < 0.999){
      continue;
    }
    double angle_dir_first_z_axis = direction_first.Angle(z_axis);
    TVector3 orthogonal_vector = direction_first.Cross(z_axis);
    if (orthogonal_vector.Unit().Mag() < 0.999){
        continue;
    }
    direction_first.Rotate(angle_dir_first_z_axis, orthogonal_vector);
    direction_second.Rotate(angle_dir_first_z_axis, orthogonal_vector);

    //Now work out the angle between the vectors in the x-z plane
    direction_first.SetY(0);
    direction_second.SetY(0);
    double dot_product = direction_first.Dot(direction_second);
    dot_product = std::min(std::max(dot_product,-1.),1.);
    double angle = acos(dot_product) * 180/3.142;

    //define +x as a +angle
    if (direction_second.X() < 0) angle*=-1;
    deflection_angles.push_back(angle);
  }

  double angle_mean = 0;
  for (size_t i_angle = 0; i_angle < deflection_angles.size(); i_angle++){
    angle_mean += deflection_angles[i_angle];
  }

  if (deflection_angles.size()>0) angle_mean/=deflection_angles.size();
  else angle_mean=-100;

  double angle_var = 0;
  for (size_t i_angle = 0; i_angle < deflection_angles.size(); i_angle++){
    angle_var = (deflection_angles[i_angle] - angle_mean)*(deflection_angles[i_angle] - angle_mean);
  }

  if (deflection_angles.size() > 1) angle_var /= (deflection_angles.size()-1);
  else angle_var = -2.;

  if (angle_var > 0.0)
      trackVars.SetWobble(sqrt(angle_var));
}
    
/////////////////////////////////////////////////////////////

void TrackVarManager::FillTrackMomentum(const art::Ptr<recob::Track> &track, 
    TrackVarManager::TrackVars &trackVars) const
{
    trkf::TrackMomentumCalculator trackMomCalc(m_minTrackLengthMCS, m_maxTrackLengthMCS);

    // Assume that it's a muon
    float byRange = trackMomCalc.GetTrackMomentum(track->Length(), 13);
    byRange = (byRange - m_intTrkMomRange) / m_gradTrkMomRange;

    float byMCS = trackMomCalc.GetMomentumMultiScatterChi2(track, true);
    byMCS = (byMCS - m_intTrkMomMCS) / m_gradTrkMomMCS;

    const float comparison = std::fabs(byRange - byMCS) / byRange;

    trackVars.SetMomentumComparison(comparison);
}

/////////////////////////////////////////////////////////////

void TrackVarManager::NormaliseTrackVars(TrackVarManager::TrackVars &trackVars) const
{
    if (trackVars.GetIsNormalised())
        throw cet::exception("ivysaur::TrackVarManager::NormaliseTrackVars: ") << "track vars are already normalised!";        

    if (trackVars.GetNTrackChildren().second)
        trackVars.SetNTrackChildren(this->NormaliseTrackVar(trackVars.GetNTrackChildren(), m_nTrackChildrenMean, m_nTrackChildrenStd));

    if (trackVars.GetNShowerChildren().second)
        trackVars.SetNShowerChildren(this->NormaliseTrackVar(trackVars.GetNShowerChildren(), m_nShowerChildrenMean, m_nShowerChildrenStd));

    if (trackVars.GetNGrandChildren().second)
        trackVars.SetNGrandChildren(this->NormaliseTrackVar(trackVars.GetNGrandChildren(), m_nGrandChildrenMean, m_nGrandChildrenStd));

    if (trackVars.GetNChildHits().second)
        trackVars.SetNChildHits(this->NormaliseTrackVar(trackVars.GetNChildHits(), m_nChildHitsMean, m_nChildHitsStd));

    if (trackVars.GetChildEnergy().second)
        trackVars.SetChildEnergy(this->NormaliseTrackVar(trackVars.GetChildEnergy(), m_childEnergyMean, m_childEnergyStd));

    if (trackVars.GetChildTrackScore().second)
        trackVars.SetChildTrackScore(this->NormaliseTrackVar(trackVars.GetChildTrackScore(), m_childTrackScoreMean, m_childTrackScoreStd));

    if (trackVars.GetTrackLength().second)
        trackVars.SetTrackLength(this->NormaliseTrackVar(trackVars.GetTrackLength(), m_trackLengthMean, m_trackLengthStd));

    if (trackVars.GetWobble().second)
        trackVars.SetWobble(this->NormaliseTrackVar(trackVars.GetWobble(), m_wobbleMean, m_wobbleStd));

    if (trackVars.GetMomentumComparison().second)
        trackVars.SetMomentumComparison(this->NormaliseTrackVar(trackVars.GetMomentumComparison(), m_momentumComparisonMean, m_momentumComparisonStd));
    
    trackVars.SetIsNormalised(true);
}
    
/////////////////////////////////////////////////////////////    
    
float TrackVarManager::NormaliseTrackVar(const std::pair<float, bool> &inputTrackVar, const float mean, const float std) const
{
    // Is still on default value...
    if (!inputTrackVar.second)
        throw cet::exception("ivysaur::TrackVarManager::NormaliseTrackVar: ") << "can't normalise a variables that wasn't set";  
    
    const float normalised((inputTrackVar.first - mean) / std);
    return normalised;
}

/////////////////////////////////////////////////////////////

}
