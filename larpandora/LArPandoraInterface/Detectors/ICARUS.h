/**
 *  @file   larpandora/LArPandoraInterface/Detectors/ICARUS.h
 *
 *  @brief  Detector interface for ICARUS
 *
 *  $Log: $
 */

#include "larpandora/LArPandoraInterface/Detectors/VintageLArTPCThreeView.h"

#include "larcore/Geometry/Geometry.h"

#include <cmath>

namespace lar_pandora {

  /**
   *  @brief  Detector interface for ICARUS
   */
  class ICARUS : public VintageLArTPCThreeView {
  public:
    geo::View_t TargetViewU(const geo::TPCID::TPCID_t tpc,
                            const geo::CryostatID::CryostatID_t cstat) const override;

    geo::View_t TargetViewV(const geo::TPCID::TPCID_t tpc,
                            const geo::CryostatID::CryostatID_t cstat) const override;

    geo::View_t TargetViewW(const geo::TPCID::TPCID_t tpc,
                            const geo::CryostatID::CryostatID_t cstat) const override;

    float WireAngleW(const geo::TPCID::TPCID_t tpc,
                     const geo::CryostatID::CryostatID_t cstat) const override;
  };

  inline geo::View_t ICARUS::TargetViewU(const geo::TPCID::TPCID_t tpc,
                                         const geo::CryostatID::CryostatID_t cstat) const
  {
    geo::TPCID const tpcID{cstat, tpc};
    return GetLArSoftGeometry().TPC(tpcID).DriftSign() == geo::DriftSign::Positive ?
             GetChannelMap().Plane(geo::PlaneID(tpcID, 1)).View() :
             GetChannelMap().Plane(geo::PlaneID(tpcID, 2)).View();
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline geo::View_t ICARUS::TargetViewV(const geo::TPCID::TPCID_t tpc,
                                         const geo::CryostatID::CryostatID_t cstat) const
  {
    geo::TPCID const tpcID{cstat, tpc};
    return GetLArSoftGeometry().TPC(tpcID).DriftSign() == geo::DriftSign::Positive ?
             GetChannelMap().Plane(geo::PlaneID(tpcID, 2)).View() :
             GetChannelMap().Plane(geo::PlaneID(tpcID, 1)).View();
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline geo::View_t ICARUS::TargetViewW(const geo::TPCID::TPCID_t tpc,
                                         const geo::CryostatID::CryostatID_t cstat) const
  {
    return GetChannelMap().Plane(geo::PlaneID(cstat, tpc, 0)).View();
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline float ICARUS::WireAngleW(const geo::TPCID::TPCID_t tpc,
                                  const geo::CryostatID::CryostatID_t cstat) const
  {
    return std::abs(
      detector_functions::WireAngle(TargetViewW(tpc, cstat), tpc, cstat, GetChannelMap()));
  }

} // namespace lar_pandora
