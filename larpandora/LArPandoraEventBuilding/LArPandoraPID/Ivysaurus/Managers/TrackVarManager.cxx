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
        m_nTrackChildren(-1),
        m_nShowerChildren(-1),
        m_nGrandChildren(-1),
        m_nChildHits(-1),
        m_childEnergy(-1.f),
        m_childTrackScore(-1.f),
        m_trackLength(-1.f),
        m_wobble(-1.f),
        m_momentumComparison(-1.f)
{
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TrackVarManager::TrackVarManager(const fhicl::ParameterSet& pset) :
    m_recoModuleLabel(pset.get<std::string>("RecoModuleLabel")),
    m_trackModuleLabel(pset.get<std::string>("TrackModuleLabel")),
    m_childSeparation(pset.get<float>("ChildSeparation")),
    m_recombFactor(pset.get<float>("RecombFactor")),
    m_calorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
    m_minTrackLengthMCS(pset.get<float>("MinTrackLengthMCS")),
    m_maxTrackLengthMCS(pset.get<float>("MaxTrackLengthMCS")),
    m_intTrkMomRange(pset.get<float>("IntTrkMomRange")),
    m_gradTrkMomRange(pset.get<float>("GradTrkMomRange")),
    m_intTrkMomMCS(pset.get<float>("IntTrkMomMCS")),
    m_gradTrkMomMCS(pset.get<float>("GradTrkMomMCS")),
    m_nTrackChildrenLimit(pset.get<float>("NTrackChildrenLimit")),
    m_nShowerChildrenLimit(pset.get<float>("NShowerChildrenLimit")),
    m_nGrandChildrenLimit(pset.get<float>("NGrandChildrenLimit")),
    m_nChildHitsLimit(pset.get<float>("NChildHitsLimit")),
    m_childEnergyLimit(pset.get<float>("ChildEnergyLimit")),
    m_childTrackScoreLimit(pset.get<float>("ChildTrackScoreLimit")),
    m_trackLengthLimit(pset.get<float>("TrackLengthLimit")),
    m_wobbleLimit(pset.get<float>("WobbleLimit")),
    m_momentumComparisonLimit(pset.get<float>("MomentumComparisonLimit"))
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

    int nTracks = 0;
    int nShowers = 0;
    int nGrandChildren = 0;

    int highestHits = -1;
    float highestHitEnergy = -1;
    float highestHitTrackScore = -1; 

    // Look at children near endpoint? 
    for (const art::Ptr<recob::PFParticle> &childPFP : childPFPs)
    {
        const std::vector<art::Ptr<recob::SpacePoint>> spacepoints = lar_pandora::PandoraPFParticleUtils::GetSpacePoints(childPFP, evt, m_recoModuleLabel);

        if (spacepoints.empty())
            continue;

        float closestSepSq = std::numeric_limits<float>::max();

        for (const art::Ptr<recob::SpacePoint> &spacepoint : spacepoints)
        {
            const TVector3 spacepointPos = TVector3(spacepoint->position().X(), spacepoint->position().Y(), spacepoint->position().Z());
            const float separationSq = (spacepointPos - trackEndpoint).Mag2();

            if (separationSq < closestSepSq)
                closestSepSq = separationSq;
        }

        if (closestSepSq > (m_childSeparation * m_childSeparation))
            continue;

        const std::vector<art::Ptr<recob::Hit>> hits = lar_pandora::PandoraPFParticleUtils::GetHits(childPFP, evt, m_recoModuleLabel);

        if (static_cast<int>(hits.size()) > highestHits)
        {
            highestHits = hits.size();
            highestHitEnergy = GetChildEnergy(evt, pfparticle);
            highestHitTrackScore = GetChildTrackScore(evt, pfparticle);
        }

        if (childPFP->PdgCode() == 13)
            ++nTracks;
        else if (childPFP->PdgCode() == 11)
            ++nShowers;

        const std::vector<art::Ptr<recob::PFParticle>> grandChildPFPs = lar_pandora::PandoraPFParticleUtils::GetChildParticles(childPFP, evt, m_recoModuleLabel);

        nGrandChildren += grandChildPFPs.size();
    }

    trackVars.SetNTrackChildren(nTracks);
    trackVars.SetNShowerChildren(nShowers);
    trackVars.SetNGrandChildren(nGrandChildren);
    trackVars.SetNChildHits(highestHits);
    trackVars.SetChildEnergy(highestHitEnergy);
    trackVars.SetChildTrackScore(highestHitTrackScore);
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

float TrackVarManager::GetChildTrackScore(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle) const
{
    const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata = lar_pandora::PandoraPFParticleUtils::GetMetadata(pfparticle, evt, m_recoModuleLabel);
    const auto metaMap = metadata->GetPropertiesMap();

    if (metaMap.find("TrackScore") == metaMap.end())
        return -1;

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
void TrackVarManager::FillWobble(const art::Ptr<recob::Track> &track, TrackVarManager::TrackVars &trackVars) const
{
  // This follows the MicroBooNE method for MCS

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
        return;

    // nTrackChildren
    if (trackVars.GetNTrackChildren() > m_nTrackChildrenLimit)
        trackVars.SetNTrackChildren(m_nTrackChildrenLimit);

    if (trackVars.GetNTrackChildren() < ((-1.0) * std::numeric_limits<float>::epsilon()))
        trackVars.SetNTrackChildren((-1.0) * m_nTrackChildrenLimit);

    trackVars.SetNTrackChildren(trackVars.GetNTrackChildren() / m_nTrackChildrenLimit);

    // nShowerChildren
    if (trackVars.GetNShowerChildren() > m_nShowerChildrenLimit)
        trackVars.SetNShowerChildren(m_nShowerChildrenLimit);

    if (trackVars.GetNShowerChildren() < ((-1.0) * std::numeric_limits<float>::epsilon()))
        trackVars.SetNShowerChildren((-1.0) * m_nShowerChildrenLimit);

    trackVars.SetNShowerChildren(trackVars.GetNShowerChildren() / m_nShowerChildrenLimit);

    // nGrandChildren
    if (trackVars.GetNGrandChildren() > m_nGrandChildrenLimit)
        trackVars.SetNGrandChildren(m_nGrandChildrenLimit);

    if (trackVars.GetNGrandChildren() < ((-1.0) * std::numeric_limits<float>::epsilon()))
        trackVars.SetNGrandChildren((-1.0) * m_nGrandChildrenLimit);

    trackVars.SetNGrandChildren(trackVars.GetNGrandChildren() / m_nGrandChildrenLimit);

    // nChildHits
    if (trackVars.GetNChildHits() > m_nChildHitsLimit)
        trackVars.SetNChildHits(m_nChildHitsLimit);

    if (trackVars.GetNChildHits() < std::numeric_limits<float>::epsilon())
        trackVars.SetNChildHits((-1.0) * m_nChildHitsLimit);

    trackVars.SetNChildHits(trackVars.GetNChildHits() / m_nChildHitsLimit);

    // childEnergy
    if (trackVars.GetChildEnergy() > m_childEnergyLimit)
        trackVars.SetChildEnergy(m_childEnergyLimit);

    if (trackVars.GetChildEnergy() < std::numeric_limits<float>::epsilon())
        trackVars.SetChildEnergy((-1.0) * m_childEnergyLimit);

    trackVars.SetChildEnergy(trackVars.GetChildEnergy() / m_childEnergyLimit);

    // childTrackScore
    if (trackVars.GetChildTrackScore() > m_childTrackScoreLimit)
        trackVars.SetChildTrackScore(m_childTrackScoreLimit);

    if (trackVars.GetChildTrackScore() < std::numeric_limits<float>::epsilon())
        trackVars.SetChildTrackScore((-1.0) * m_childTrackScoreLimit);

    trackVars.SetChildTrackScore(trackVars.GetChildTrackScore() / m_childTrackScoreLimit);

    // trackLength
    if (trackVars.GetTrackLength() > m_trackLengthLimit)
        trackVars.SetTrackLength(m_trackLengthLimit);

    if (trackVars.GetTrackLength() < std::numeric_limits<float>::epsilon())
        trackVars.SetTrackLength((-1.0) * m_trackLengthLimit);

    trackVars.SetTrackLength(trackVars.GetTrackLength() / m_trackLengthLimit);

    // wobble
    if (trackVars.GetWobble() > m_wobbleLimit)
        trackVars.SetWobble(m_wobbleLimit);

    if (trackVars.GetWobble() < std::numeric_limits<float>::epsilon())
        trackVars.SetWobble((-1.0) * m_wobbleLimit);

    trackVars.SetWobble(trackVars.GetWobble() / m_wobbleLimit);

    // momentumComparison
    if (trackVars.GetMomentumComparison() > m_momentumComparisonLimit)
        trackVars.SetMomentumComparison(m_momentumComparisonLimit);

    if (trackVars.GetMomentumComparison() < std::numeric_limits<float>::epsilon())
        trackVars.SetMomentumComparison((-1.0) * m_momentumComparisonLimit);

    trackVars.SetMomentumComparison(trackVars.GetMomentumComparison() / m_momentumComparisonLimit);

    // Set normalised
    trackVars.SetIsNormalised(true);
}

/////////////////////////////////////////////////////////////

void TrackVarManager::Reset(TrackVarManager::TrackVars &trackVars) const
{
    trackVars.SetNTrackChildren(-1);
    trackVars.SetNShowerChildren(-1);
    trackVars.SetNGrandChildren(-1);
    trackVars.SetNChildHits(-1);
    trackVars.SetChildEnergy(-1);
    trackVars.SetChildTrackScore(-1);
    trackVars.SetTrackLength(-1.f);
    trackVars.SetWobble(-1.f);
    trackVars.SetMomentumComparison(-1.f);
}

/////////////////////////////////////////////////////////////

}
