/**
 *  @file   larpandora/LArPandoraInterface/Detectors/ProtoDUNEDualPhase.h
 *
 *  @brief  Detector interface for ProtoDUNE dual phase 
 *
 *  $Log: $
 */

#include "larpandora/LArPandoraInterface/Detectors/VintageLArTPCThreeView.h"
#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

#include "larcore/Geometry/Geometry.h"

namespace lar_pandora {

  /**
     *  @brief  Detector interface for ProtoDUNE dual phase
     */
  class ProtoDUNEDualPhase : public VintageLArTPCThreeView {
  public:
    geo::View_t TargetViewU(const geo::TPCID::TPCID_t tpc,
                            const geo::CryostatID::CryostatID_t cstat) const override;

    geo::View_t TargetViewV(const geo::TPCID::TPCID_t tpc,
                            const geo::CryostatID::CryostatID_t cstat) const override;

    geo::View_t TargetViewW(const geo::TPCID::TPCID_t tpc,
                            const geo::CryostatID::CryostatID_t cstat) const override;

    float WirePitchW() const override;

    float WireAngleW(const geo::TPCID::TPCID_t tpc,
                     const geo::CryostatID::CryostatID_t cstat) const override;

    bool CheckDetectorGapSize(const geo::Vector_t& gaps,
                              const geo::Vector_t& deltas,
                              const float maxDisplacement) const override;

    LArDetectorGap CreateDetectorGap(const geo::Point_t& point1,
                                     const geo::Point_t& point2,
                                     const geo::Vector_t& widths) const override;

    void LoadDaughterDetectorGaps(const LArDriftVolume& driftVolume,
                                  const float maxDisplacement,
                                  LArDetectorGapList& listOfGaps) const override;

    PandoraApi::Geometry::LineGap::Parameters CreateLineGapParametersFromDetectorGaps(
      const LArDetectorGap& gap) const override;

    PandoraApi::Geometry::LineGap::Parameters CreateLineGapParametersFromReadoutGaps(
      const geo::View_t view,
      const geo::TPCID::TPCID_t tpc,
      const geo::CryostatID::CryostatID_t cstat,
      const double firstXYZ[3],
      const double lastXYZ[3],
      const float halfWirePitch,
      const float xFirst,
      const float xLast,
      const pandora::Pandora* pPandora) const override;
  };

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline geo::View_t
  ProtoDUNEDualPhase::TargetViewU(const geo::TPCID::TPCID_t tpc,
                                  const geo::CryostatID::CryostatID_t cstat) const
  {
    return this->GetLArSoftGeometry()->View(geo::PlaneID(cstat, tpc, 1));
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline geo::View_t
  ProtoDUNEDualPhase::TargetViewV(const geo::TPCID::TPCID_t tpc,
                                  const geo::CryostatID::CryostatID_t cstat) const
  {
    return this->GetLArSoftGeometry()->View(geo::PlaneID(cstat, tpc, 0));
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline geo::View_t
  ProtoDUNEDualPhase::TargetViewW(const geo::TPCID::TPCID_t tpc,
                                  const geo::CryostatID::CryostatID_t cstat) const
  {
    return geo::kUnknown;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline float
  ProtoDUNEDualPhase::WirePitchW() const
  {
    return 0.5f * (this->WirePitchU() + this->WirePitchV());
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline float
  ProtoDUNEDualPhase::WireAngleW(const geo::TPCID::TPCID_t tpc,
                                 const geo::CryostatID::CryostatID_t cstat) const
  {
    return std::numeric_limits<float>::epsilon();
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline bool
  ProtoDUNEDualPhase::CheckDetectorGapSize(const geo::Vector_t& gaps,
                                           const geo::Vector_t& deltas,
                                           const float maxDisplacement) const
  {
    return (std::fabs(gaps.Y()) > maxDisplacement || std::fabs(gaps.Z()) > maxDisplacement);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline LArDetectorGap
  ProtoDUNEDualPhase::CreateDetectorGap(const geo::Point_t& point1,
                                        const geo::Point_t& point2,
                                        const geo::Vector_t& widths) const
  {
    return LArDetectorGap(point1.X(),
                          point1.Y() + widths.Y(),
                          point1.Z() + widths.Z(),
                          point2.X(),
                          point2.Y() - widths.Y(),
                          point2.Z() - widths.Z());
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline PandoraApi::Geometry::LineGap::Parameters
  ProtoDUNEDualPhase::CreateLineGapParametersFromDetectorGaps(const LArDetectorGap& gap) const
  {
    const bool isGapInU((
      std::fabs(gap.GetY2() - gap.GetY1()) >
      gap
        .GetMaxGapSize())); //Could have chosen Z here, resulting in switching Y<->Z and U<->V in the block below
    PandoraApi::Geometry::LineGap::Parameters parameters;
    parameters.m_lineGapType =
      (isGapInU ?
         pandora::TPC_WIRE_GAP_VIEW_U :
         pandora::
           TPC_WIRE_GAP_VIEW_V); //If gapSizeY is too large then the gap is in Z, therefore should be in kU (i.e. kZ)
    parameters.m_lineStartX = gap.GetX2();
    parameters.m_lineEndX = gap.GetX1();
    parameters.m_lineStartZ = (isGapInU ? gap.GetZ1() : gap.GetY1());
    parameters.m_lineEndZ = (isGapInU ? gap.GetZ2() : gap.GetY2());

    return parameters;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline PandoraApi::Geometry::LineGap::Parameters
  ProtoDUNEDualPhase::CreateLineGapParametersFromReadoutGaps(
    const geo::View_t view,
    const geo::TPCID::TPCID_t tpc,
    const geo::CryostatID::CryostatID_t cstat,
    const double firstXYZ[3],
    const double lastXYZ[3],
    const float halfWirePitch,
    const float xFirst,
    const float xLast,
    const pandora::Pandora* pPandora) const
  {
    float first(0.f), last(0.f);
    pandora::LineGapType gapType(pandora::TPC_DRIFT_GAP);
    if (view == this->TargetViewU(tpc, cstat)) {
      first = firstXYZ[2];
      last = lastXYZ[2];
      gapType = pandora::TPC_WIRE_GAP_VIEW_U;
    }
    else if (view == this->TargetViewV(tpc, cstat)) {
      first = firstXYZ[1];
      last = lastXYZ[1];
      gapType = pandora::TPC_WIRE_GAP_VIEW_V;
    }
    return detector_functions::CreateReadoutGapParameters(
      first, last, xFirst, xLast, halfWirePitch, gapType);
  }

} // namespace lar_pandora
