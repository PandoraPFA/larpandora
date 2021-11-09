/**
 *  @file   larpandora/LArPandoraInterface/Detectors/ICARUS.h
 *
 *  @brief  Detector interface for ICARUS
 *
 *  $Log: $
 */

#include "larpandora/LArPandoraInterface/Detectors/VintageLArTPCThreeView.h"
#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

#include "larcore/Geometry/Geometry.h"

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

  inline geo::View_t
  ICARUS::TargetViewU(const geo::TPCID::TPCID_t tpc,
                      const geo::CryostatID::CryostatID_t cstat) const
  {
    return (this->GetLArSoftGeometry()->TPC(tpc, cstat).DriftDirection() == geo::kPosX ?
              this->GetLArSoftGeometry()->View(geo::PlaneID(cstat, tpc, 1)) :
              this->GetLArSoftGeometry()->View(geo::PlaneID(cstat, tpc, 2)));
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline geo::View_t
  ICARUS::TargetViewV(const geo::TPCID::TPCID_t tpc,
                      const geo::CryostatID::CryostatID_t cstat) const
  {
    return (this->GetLArSoftGeometry()->TPC(tpc, cstat).DriftDirection() == geo::kPosX ?
              this->GetLArSoftGeometry()->View(geo::PlaneID(cstat, tpc, 2)) :
              this->GetLArSoftGeometry()->View(geo::PlaneID(cstat, tpc, 1)));
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline geo::View_t
  ICARUS::TargetViewW(const geo::TPCID::TPCID_t tpc,
                      const geo::CryostatID::CryostatID_t cstat) const
  {
    return this->GetLArSoftGeometry()->View(geo::PlaneID(cstat, tpc, 0));
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline float
  ICARUS::WireAngleW(const geo::TPCID::TPCID_t tpc, const geo::CryostatID::CryostatID_t cstat) const
  {
    return std::abs(detector_functions::WireAngle(
      this->TargetViewW(tpc, cstat), tpc, cstat, this->GetLArSoftGeometry()));
  }

} // namespace lar_pandora
