#ifndef IVYSAURUSUTILS_H_SEEN
#define IVYSAURUSUTILS_H_SEEN

///////////////////////////////////////////////
// IvysaurusUtils.h
//
// Isobel Mawby i.mawby1@lancaster.ac.uk
///////////////////////////////////////////////

// Framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larsim/Utils/TruthMatchUtils.h"

// c++
#include <vector>
#include <map>

namespace IvysaurusUtils
{
    enum PandoraView {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W};

    // Function to project a 3D coordinate into a specified Pandora 2D view
    const TVector3 ProjectIntoPandoraView(const TVector3 &inputPosition3D, const PandoraView pandoraView);

    // Function to obtain the Pandora view of a LArSoft hit
    const PandoraView GetPandora2DView(const art::Ptr<recob::Hit> &hit);

    void ObtainPandoraHitPositionAndWidth(const art::Event &evt, const art::Ptr<recob::Hit> hit, 
        const PandoraView hitType, TVector3 &pandoraHitPosition, float &width);

    // Function to obtain the Pandora U coordinate from LArSoft Y/Z coordinates
    float YZToU(const float yCoord, const float zCoord);

    // Function to obtain the Pandora V coordinate from LArSoft Y/Z coordinates
    float YZToV(const float yCoord, const float zCoord);

    // Function to obtain the Pandora W coordinate from LArSoft Y/Z coordinates
    float YZToW(const float yCoord, const float zCoord);

    double CompletenessFromTrueParticleID(detinfo::DetectorClocksData const& clockData, const std::vector<art::Ptr<recob::Hit>> &selectedHits, 
        const std::vector<art::Ptr<recob::Hit>> &eventHits, const int trackID);

    double HitPurityFromTrueParticleID(detinfo::DetectorClocksData const& clockData, const std::vector<art::Ptr<recob::Hit>> &selectedHits,
        const int trackID);

    float IntegrateGaussian(const float limitA, const float limitB, const float mean, const float std, const float stepSize);

    float TrapeziumRule(const float lowerLimit, const float upperLimit, const float mean, const float std);

    bool GetInitialDirection(const art::Event &evt, const TVector3 &pfpVertex, const std::vector<art::Ptr<recob::SpacePoint>> &spacepoints, 
        const std::string &recoModuleLabel, TVector3 &direction);
}

#endif
