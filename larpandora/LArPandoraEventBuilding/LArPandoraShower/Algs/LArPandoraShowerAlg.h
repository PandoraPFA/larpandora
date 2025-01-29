#ifndef LArPandoraShowerAlg_hxx
#define LArPandoraShowerAlg_hxx

namespace reco::shower {
  class ShowerElementHolder;
}

namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

#include "larcore/Geometry/WireReadout.h"

namespace recob {
  class Hit;
  class PFParticle;
  class SpacePoint;
}

namespace spacecharge {
  class SpaceCharge;
}

#include "art/Framework/Principal/fwd.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"

//C++ Includes
#include <string>
#include <vector>

namespace shower {
  class LArPandoraShowerAlg;
}

class shower::LArPandoraShowerAlg {
public:
  explicit LArPandoraShowerAlg(const fhicl::ParameterSet& pset);

  void OrderShowerHits(detinfo::DetectorPropertiesData const& detProp,
                       std::vector<art::Ptr<recob::Hit>>& hits,
                       geo::Point_t const& ShowerPosition,
                       geo::Vector_t const& ShowerDirection) const;

  void OrderShowerSpacePointsPerpendicular(std::vector<art::Ptr<recob::SpacePoint>>& showersps,
                                           geo::Point_t const& vertex,
                                           geo::Vector_t const& direction) const;

  void OrderShowerSpacePoints(std::vector<art::Ptr<recob::SpacePoint>>& showersps,
                              geo::Point_t const& vertex,
                              geo::Vector_t const& direction) const;

  void OrderShowerSpacePoints(std::vector<art::Ptr<recob::SpacePoint>>& showersps,
                              geo::Point_t const& vertex) const;

  geo::Point_t ShowerCentre(std::vector<art::Ptr<recob::SpacePoint>> const& showersps) const;

  geo::Point_t ShowerCentre(detinfo::DetectorClocksData const& clockData,
                            detinfo::DetectorPropertiesData const& detProp,
                            std::vector<art::Ptr<recob::SpacePoint>> const& showersps,
                            art::FindManyP<recob::Hit> const& fmh,
                            float& totalCharge) const;

  geo::Point_t ShowerCentre(detinfo::DetectorClocksData const& clockData,
                            detinfo::DetectorPropertiesData const& detProp,
                            std::vector<art::Ptr<recob::SpacePoint>> const& showerspcs,
                            art::FindManyP<recob::Hit> const& fmh) const;

  double DistanceBetweenSpacePoints(art::Ptr<recob::SpacePoint> const& sp_a,
                                    art::Ptr<recob::SpacePoint> const& sp_b) const;

  double SpacePointCharge(art::Ptr<recob::SpacePoint> const& sp,
                          art::FindManyP<recob::Hit> const& fmh) const;

  double SpacePointTime(art::Ptr<recob::SpacePoint> const& sp,
                        art::FindManyP<recob::Hit> const& fmh) const;

  TVector2 HitCoordinates(detinfo::DetectorPropertiesData const& detProp,
                          art::Ptr<recob::Hit> const& hit) const;

  double SpacePointProjection(art::Ptr<recob::SpacePoint> const& sp,
                              geo::Point_t const& vertex,
                              geo::Vector_t const& direction) const;

  double SpacePointPerpendicular(art::Ptr<recob::SpacePoint> const& sp,
                                 geo::Point_t const& vertex,
                                 geo::Vector_t const& direction) const;

  double SpacePointPerpendicular(art::Ptr<recob::SpacePoint> const& sp,
                                 geo::Point_t const& vertex,
                                 geo::Vector_t const& direction,
                                 double proj) const;

  double RMSShowerGradient(std::vector<art::Ptr<recob::SpacePoint>>& sps,
                           const geo::Point_t& ShowerCentre,
                           const geo::Vector_t& Direction,
                           const unsigned int nSegments) const;

  double CalculateRMS(const std::vector<float>& perps) const;

  // The SCE service requires thing in geo::Point/Vector form, so overload and be nice
  double SCECorrectPitch(double const& pitch,
                         geo::Point_t const& pos,
                         geo::Vector_t const& dir,
                         unsigned int const& TPC) const;

  double SCECorrectEField(double const& EField, geo::Point_t const& pos) const;

  std::map<art::Ptr<recob::Hit>, std::vector<art::Ptr<recob::Hit>>> OrganizeHits(
    const std::vector<art::Ptr<recob::Hit>>& hits) const;

  void DebugEVD(art::Ptr<recob::PFParticle> const& pfparticle,
                art::Event const& Event,
                const reco::shower::ShowerElementHolder& ShowerEleHolder,
                std::string const& evd_disp_name_append = "") const;

private:
  bool fUseCollectionOnly;
  art::InputTag fPFParticleLabel;
  bool fSCEXFlip; // If a (legacy) flip is needed in x componant of spatial SCE correction

  spacecharge::SpaceCharge const* fSCE;
  geo::WireReadoutGeom const* fChannelMap = &art::ServiceHandle<geo::WireReadout>()->Get();
  art::ServiceHandle<art::TFileService> tfs;

  const std::string fInitialTrackInputLabel;
  const std::string fShowerStartPositionInputLabel;
  const std::string fShowerDirectionInputLabel;
  const std::string fInitialTrackSpacePointsInputLabel;
};

#endif
