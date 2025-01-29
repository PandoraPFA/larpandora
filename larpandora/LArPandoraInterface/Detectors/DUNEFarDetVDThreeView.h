/**
 *  @file   larpandora/LArPandoraInterface/Detectors/DUNEFarDetVDThreeView.h
 *
 *  @brief  Detector interface for DUNE's vertical drift, 3 view far detector
 *
 *  $Log: $
 */

#include "larpandora/LArPandoraInterface/Detectors/VintageLArTPCThreeView.h"

#include "larcore/Geometry/Geometry.h"

namespace lar_pandora {

  /**
     *  @brief  Detector interface DUNE's vertical drift far detector
     */
  class DUNEFarDetVDThreeView : public VintageLArTPCThreeView {
  public:
    geo::View_t TargetViewU(const geo::TPCID::TPCID_t tpc,
                            const geo::CryostatID::CryostatID_t cstat) const override;

    geo::View_t TargetViewV(const geo::TPCID::TPCID_t tpc,
                            const geo::CryostatID::CryostatID_t cstat) const override;

    geo::View_t TargetViewW(const geo::TPCID::TPCID_t tpc,
                            const geo::CryostatID::CryostatID_t cstat) const override;
  };

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline geo::View_t DUNEFarDetVDThreeView::TargetViewU(
    const geo::TPCID::TPCID_t tpc,
    const geo::CryostatID::CryostatID_t cstat) const
  {
    return this->GetChannelMap().Plane(geo::PlaneID(cstat, tpc, 0)).View();
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline geo::View_t DUNEFarDetVDThreeView::TargetViewV(
    const geo::TPCID::TPCID_t tpc,
    const geo::CryostatID::CryostatID_t cstat) const
  {
    return this->GetChannelMap().Plane(geo::PlaneID(cstat, tpc, 1)).View();
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline geo::View_t DUNEFarDetVDThreeView::TargetViewW(
    const geo::TPCID::TPCID_t tpc,
    const geo::CryostatID::CryostatID_t cstat) const
  {
    return this->GetChannelMap().Plane(geo::PlaneID(cstat, tpc, 2)).View();
  }

} // namespace lar_pandora
