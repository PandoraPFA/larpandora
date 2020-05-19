#ifndef LArPandoraShowerCheatingAlg_hxx
#define LArPandoraShowerCheatingAlg_hxx

//Framework Includes
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft Includes
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Algs/ShowerElementHolder.hh"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Algs/LArPandoraShowerAlg.h"

//C++ Includes
#include <iostream>
#include <vector>
#include <map>

//Root Includes
#include "TVector3.h"
#include "TMath.h"
#include "TVector.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TString.h"

namespace shower {
  class LArPandoraShowerCheatingAlg;
}

class shower::LArPandoraShowerCheatingAlg {
  public:
    LArPandoraShowerCheatingAlg(const fhicl::ParameterSet& pset);

    std::map<int,const simb::MCParticle*> GetTrueParticleMap() const;
    std::map<int,std::vector<int> > GetTrueChain(std::map<int,const simb::MCParticle*>& trueParticles) const;
    void CheatDebugEVD(const simb::MCParticle* trueParticle, art::Event const& Event,
        reco::shower::ShowerElementHolder& ShowerEleHolder,
        const art::Ptr<recob::PFParticle>& pfparticle) const;

    int TrueParticleID(const art::Ptr<recob::Hit>& hit) const;

    std::pair<int,double> TrueParticleIDFromTrueChain(std::map<int,std::vector<int> > const& ShowersMothers,
                                                      std::vector<art::Ptr<recob::Hit> > const& hits, int planeid) const;

  private:

    shower::LArPandoraShowerAlg fLArPandoraShowerAlg;

    art::InputTag                                       fHitModuleLabel;
    art::InputTag                                       fPFParticleLabel;
    art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
    art::ServiceHandle<art::TFileService>   tfs;


    std::string fShowerStartPositionInputLabel;
    std::string fShowerDirectionInputLabel;
    std::string fInitialTrackSpacePointsInputLabel;

};
#endif
