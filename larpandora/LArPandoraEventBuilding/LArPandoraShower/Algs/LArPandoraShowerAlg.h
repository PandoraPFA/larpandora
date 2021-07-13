#ifndef LArPandoraShowerAlg_hxx
#define LArPandoraShowerAlg_hxx

//Framework Includes
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//LArSoft Includes
#include "art/Framework/Principal/Event.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Algs/ShowerElementHolder.hh"

//C++ Includes
#include <algorithm>
#include <iostream>
#include <map>
#include <vector>

//Root Includes
#include "TCanvas.h"
#include "TH3F.h"
#include "TMath.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include "TVector.h"
#include "TVector3.h"

namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

namespace shower {
  class LArPandoraShowerAlg;
}

class shower::LArPandoraShowerAlg {
public:
  explicit LArPandoraShowerAlg(const fhicl::ParameterSet& pset);

  void OrderShowerHits(detinfo::DetectorPropertiesData const& detProp,
                       std::vector<art::Ptr<recob::Hit>>& hits,
                       TVector3 const& ShowerDirection,
                       TVector3 const& ShowerPosition) const;

  void OrderShowerSpacePointsPerpendicular(std::vector<art::Ptr<recob::SpacePoint>>& showersps,
                                           TVector3 const& vertex,
                                           TVector3 const& direction) const;

  void OrderShowerSpacePoints(std::vector<art::Ptr<recob::SpacePoint>>& showersps,
                              TVector3 const& vertex,
                              TVector3 const& direction) const;

  void OrderShowerSpacePoints(std::vector<art::Ptr<recob::SpacePoint>>& showersps,
                              TVector3 const& vertex) const;

  TVector3 ShowerCentre(std::vector<art::Ptr<recob::SpacePoint>> const& showersps) const;

  TVector3 ShowerCentre(detinfo::DetectorClocksData const& clockData,
                        detinfo::DetectorPropertiesData const& detProp,
                        std::vector<art::Ptr<recob::SpacePoint>> const& showersps,
                        art::FindManyP<recob::Hit> const& fmh,
                        float& totalCharge) const;

  TVector3 ShowerCentre(detinfo::DetectorClocksData const& clockData,
                        detinfo::DetectorPropertiesData const& detProp,
                        std::vector<art::Ptr<recob::SpacePoint>> const& showerspcs,
                        art::FindManyP<recob::Hit> const& fmh) const;

  TVector3 SpacePointPosition(art::Ptr<recob::SpacePoint> const& sp) const;

  double DistanceBetweenSpacePoints(art::Ptr<recob::SpacePoint> const& sp_a,
                                    art::Ptr<recob::SpacePoint> const& sp_b) const;

  double SpacePointCharge(art::Ptr<recob::SpacePoint> const& sp,
                          art::FindManyP<recob::Hit> const& fmh) const;

  double SpacePointTime(art::Ptr<recob::SpacePoint> const& sp,
                        art::FindManyP<recob::Hit> const& fmh) const;

  TVector2 HitCoordinates(detinfo::DetectorPropertiesData const& detProp,
                          art::Ptr<recob::Hit> const& hit) const;

  double SpacePointProjection(art::Ptr<recob::SpacePoint> const& sp,
                              TVector3 const& vertex,
                              TVector3 const& direction) const;

  double SpacePointPerpendicular(art::Ptr<recob::SpacePoint> const& sp,
                                 TVector3 const& vertex,
                                 TVector3 const& direction) const;

  double SpacePointPerpendicular(art::Ptr<recob::SpacePoint> const& sp,
                                 TVector3 const& vertex,
                                 TVector3 const& direction,
                                 double proj) const;

  double RMSShowerGradient(std::vector<art::Ptr<recob::SpacePoint>>& sps,
                           const TVector3& ShowerCentre,
                           const TVector3& Direction,
                           const unsigned int nSegments) const;

  double CalculateRMS(const std::vector<float>& perps) const;

  // The SCE service requires thing in geo::Point/Vector form, so overload and be nice
  double SCECorrectPitch(double const& pitch,
                         TVector3 const& pos,
                         TVector3 const& dir,
                         unsigned int const& TPC) const;
  double SCECorrectPitch(double const& pitch,
                         geo::Point_t const& pos,
                         geo::Vector_t const& dir,
                         unsigned int const& TPC) const;

  double SCECorrectEField(double const& EField, TVector3 const& pos) const;
  double SCECorrectEField(double const& EField, geo::Point_t const& pos) const;

  void DebugEVD(art::Ptr<recob::PFParticle> const& pfparticle,
                art::Event const& Event,
                const reco::shower::ShowerElementHolder& ShowerEleHolder,
                std::string const& evd_disp_name_append = "") const;

private:
  bool fUseCollectionOnly;
  art::InputTag fPFParticleLabel;
  bool fSCEXFlip; // If a (legacy) flip is needed in x componant of spatial SCE correction

  spacecharge::SpaceCharge const* fSCE;
  art::ServiceHandle<geo::Geometry const> fGeom;
  art::ServiceHandle<art::TFileService> tfs;

  const std::string fInitialTrackInputLabel;
  const std::string fShowerStartPositionInputLabel;
  const std::string fShowerDirectionInputLabel;
  const std::string fInitialTrackSpacePointsInputLabel;
};

#endif
