/**
 *  @file   larpandora/LArPandoraInterface/Detectors/VintageLArTPCThreeView.h
 *
 *  @brief  Detector interface for an older-style 3view, horizontal drift, single-phase LArTPC (e.g. MicroBooNE, DUNEFD-HD, PDSP etc)
 *
 *  $Log: $
 */

#ifndef LAR_PANDORA_DETECTOR_VINTAGE_LAR_TPC_THREE_VIEW_H
#define LAR_PANDORA_DETECTOR_VINTAGE_LAR_TPC_THREE_VIEW_H 1

#include "larpandora/LArPandoraInterface/Detectors/LArPandoraDetectorType.h"
#include "larpandora/LArPandoraInterface/LArPandoraGeometryComponents.h"

#include "Api/PandoraApi.h"
#include "Managers/PluginManager.h"
#include "Pandora/Pandora.h"
#include "Plugins/LArTransformationPlugin.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"

namespace lar_pandora {

  /**
   *  @brief  Detector interface for a 3view, horizontal drift, single-phase LArTPC
   */
  class VintageLArTPCThreeView : public LArPandoraDetectorType {
  public:
    geo::View_t TargetViewU(const geo::TPCID::TPCID_t tpc,
                            const geo::CryostatID::CryostatID_t cstat) const override;

    geo::View_t TargetViewV(const geo::TPCID::TPCID_t tpc,
                            const geo::CryostatID::CryostatID_t cstat) const override;

    geo::View_t TargetViewW(const geo::TPCID::TPCID_t tpc,
                            const geo::CryostatID::CryostatID_t cstat) const override;

    float WirePitchU() const override;

    float WirePitchV() const override;

    float WirePitchW() const override;

    float WireAngleU(const geo::TPCID::TPCID_t tpc,
                     const geo::CryostatID::CryostatID_t cstat) const override;

    float WireAngleV(const geo::TPCID::TPCID_t tpc,
                     const geo::CryostatID::CryostatID_t cstat) const override;

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
      const geo::Point_t& firstXYZ,
      const geo::Point_t& lastXYZ,
      const float halfWirePitch,
      const float xFirst,
      const float xLast,
      const pandora::Pandora* pPandora) const override;

    /**
     *  @brief  Loan the LArSoft geometry handle owned by this class
     *
     *  @result The LArSoft geometry handle
     */
    const geo::GeometryCore& GetLArSoftGeometry() const;
    const geo::WireReadoutGeom& GetChannelMap() const;

  private:
    art::ServiceHandle<geo::Geometry> m_LArSoftGeometry; ///< the LArSoft geometry handle
    geo::WireReadoutGeom const* m_wireReadoutGeom = &art::ServiceHandle<geo::WireReadout>()->Get();
  };

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline geo::View_t VintageLArTPCThreeView::TargetViewU(
    const geo::TPCID::TPCID_t tpc,
    const geo::CryostatID::CryostatID_t cstat) const
  {
    geo::TPCID const tpcID{cstat, tpc};
    auto const [_, sign] = m_LArSoftGeometry->TPC(tpcID).DriftAxisWithSign();
    return (sign == geo::DriftSign::Positive ?
              m_wireReadoutGeom->Plane(geo::PlaneID(tpcID, 1)).View() :
              m_wireReadoutGeom->Plane(geo::PlaneID(tpcID, 0)).View());
  }

  inline geo::View_t VintageLArTPCThreeView::TargetViewV(
    const geo::TPCID::TPCID_t tpc,
    const geo::CryostatID::CryostatID_t cstat) const
  {
    geo::TPCID const tpcID{cstat, tpc};
    auto const [_, sign] = m_LArSoftGeometry->TPC(tpcID).DriftAxisWithSign();
    return (sign == geo::DriftSign::Positive ?
              m_wireReadoutGeom->Plane(geo::PlaneID(tpcID, 0)).View() :
              m_wireReadoutGeom->Plane(geo::PlaneID(tpcID, 1)).View());
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline geo::View_t VintageLArTPCThreeView::TargetViewW(
    const geo::TPCID::TPCID_t tpc,
    const geo::CryostatID::CryostatID_t cstat) const
  {
    return m_wireReadoutGeom->Plane({cstat, tpc, 2}).View();
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline float VintageLArTPCThreeView::WirePitchU() const
  {
    return m_wireReadoutGeom->Plane({0, 0}, TargetViewU(0, 0)).WirePitch();
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline float VintageLArTPCThreeView::WirePitchV() const
  {
    return m_wireReadoutGeom->Plane({0, 0}, TargetViewV(0, 0)).WirePitch();
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline float VintageLArTPCThreeView::WirePitchW() const
  {
    return m_wireReadoutGeom->Plane({0, 0}, TargetViewW(0, 0)).WirePitch();
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline float VintageLArTPCThreeView::WireAngleU(const geo::TPCID::TPCID_t tpc,
                                                  const geo::CryostatID::CryostatID_t cstat) const
  {
    return detector_functions::WireAngle(TargetViewU(tpc, cstat), tpc, cstat, *m_wireReadoutGeom);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline float VintageLArTPCThreeView::WireAngleV(const geo::TPCID::TPCID_t tpc,
                                                  const geo::CryostatID::CryostatID_t cstat) const
  {
    return detector_functions::WireAngle(TargetViewV(tpc, cstat), tpc, cstat, *m_wireReadoutGeom);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline float VintageLArTPCThreeView::WireAngleW(const geo::TPCID::TPCID_t tpc,
                                                  const geo::CryostatID::CryostatID_t cstat) const
  {
    return detector_functions::WireAngle(TargetViewW(tpc, cstat), tpc, cstat, *m_wireReadoutGeom);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline bool VintageLArTPCThreeView::CheckDetectorGapSize(const geo::Vector_t& gaps,
                                                           const geo::Vector_t& deltas,
                                                           const float maxDisplacement) const
  {
    if (gaps.X() < 0.f || gaps.X() > maxDisplacement || deltas.Y() > maxDisplacement ||
        deltas.Z() > maxDisplacement)
      return false;
    return true;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline LArDetectorGap VintageLArTPCThreeView::CreateDetectorGap(const geo::Point_t& point1,
                                                                  const geo::Point_t& point2,
                                                                  const geo::Vector_t& widths) const
  {
    return LArDetectorGap(point1.X(), point1.Y(), point1.Z(), point2.X(), point2.Y(), point2.Z());
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline void VintageLArTPCThreeView::LoadDaughterDetectorGaps(const LArDriftVolume& driftVolume,
                                                               const float maxDisplacement,
                                                               LArDetectorGapList& listOfGaps) const
  {
    return;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline PandoraApi::Geometry::LineGap::Parameters
  VintageLArTPCThreeView::CreateLineGapParametersFromDetectorGaps(const LArDetectorGap& gap) const
  {
    return detector_functions::CreateDriftGapParameters(gap);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline PandoraApi::Geometry::LineGap::Parameters
  VintageLArTPCThreeView::CreateLineGapParametersFromReadoutGaps(
    const geo::View_t view,
    const geo::TPCID::TPCID_t tpc,
    const geo::CryostatID::CryostatID_t cstat,
    const geo::Point_t& firstXYZ,
    const geo::Point_t& lastXYZ,
    const float halfWirePitch,
    const float xFirst,
    const float xLast,
    const pandora::Pandora* pPandora) const
  {
    float first(0.f), last(0.f);
    pandora::LineGapType gapType(pandora::TPC_DRIFT_GAP);
    if (view == TargetViewW(tpc, cstat)) {
      first = firstXYZ.Z();
      last = lastXYZ.Z();
      gapType = pandora::TPC_WIRE_GAP_VIEW_W;
    }
    else if (view == TargetViewU(tpc, cstat)) {
      first =
        pPandora->GetPlugins()->GetLArTransformationPlugin()->YZtoU(firstXYZ.Y(), firstXYZ.Z());
      last = pPandora->GetPlugins()->GetLArTransformationPlugin()->YZtoU(lastXYZ.Y(), lastXYZ.Z());
      gapType = pandora::TPC_WIRE_GAP_VIEW_U;
    }
    else if (view == TargetViewV(tpc, cstat)) {
      first =
        pPandora->GetPlugins()->GetLArTransformationPlugin()->YZtoV(firstXYZ.Y(), firstXYZ.Z());
      last = pPandora->GetPlugins()->GetLArTransformationPlugin()->YZtoV(lastXYZ.Y(), lastXYZ.Z());
      gapType = pandora::TPC_WIRE_GAP_VIEW_V;
    }
    return detector_functions::CreateReadoutGapParameters(
      first, last, xFirst, xLast, halfWirePitch, gapType);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline const geo::GeometryCore& VintageLArTPCThreeView::GetLArSoftGeometry() const
  {
    return *m_LArSoftGeometry;
  }

  inline const geo::WireReadoutGeom& VintageLArTPCThreeView::GetChannelMap() const
  {
    return *m_wireReadoutGeom;
  }

} // namespace lar_pandora
#endif // #ifndef LAR_PANDORA_DETECTOR_VINTAGE_LAR_TPC_THREE_VIEW_H
