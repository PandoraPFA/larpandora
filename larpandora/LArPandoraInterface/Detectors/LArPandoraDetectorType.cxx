/**
 *  @file   larpandora/LArPandoraInterface/Detectors/LArPandoraDetectorType.cxx
 *
 *  @brief  Implementation of the interface for handling detector-specific details, as well as some helper functions
 *
 *  $Log: $
 */

#include "larpandora/LArPandoraInterface/Detectors/LArPandoraDetectorType.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "larpandora/LArPandoraInterface/LArPandoraGeometryComponents.h"

#include <limits>

namespace lar_pandora {

  float detector_functions::WireAngle(const geo::View_t view,
                                      const geo::TPCID::TPCID_t tpc,
                                      const geo::CryostatID::CryostatID_t cstat,
                                      const geo::WireReadoutGeom& wireReadoutGeom)
  {
    return 0.5f * M_PI - wireReadoutGeom.WireAngleToVertical(view, geo::TPCID{cstat, tpc});
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  PandoraApi::Geometry::LineGap::Parameters detector_functions::CreateDriftGapParameters(
    const LArDetectorGap& gap)
  {
    PandoraApi::Geometry::LineGap::Parameters parameters;
    parameters.m_lineGapType = pandora::TPC_DRIFT_GAP;
    parameters.m_lineStartX = gap.GetX1();
    parameters.m_lineEndX = gap.GetX2();
    parameters.m_lineStartZ = -std::numeric_limits<float>::max();
    parameters.m_lineEndZ = std::numeric_limits<float>::max();

    return parameters;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  PandoraApi::Geometry::LineGap::Parameters detector_functions::CreateReadoutGapParameters(
    const float firstPoint,
    const float lastPoint,
    const float xFirst,
    const float xLast,
    const float halfWirePitch,
    const pandora::LineGapType gapType)
  {
    PandoraApi::Geometry::LineGap::Parameters parameters;
    parameters.m_lineGapType = gapType;
    parameters.m_lineStartX = xFirst;
    parameters.m_lineEndX = xLast;
    parameters.m_lineStartZ = std::min(firstPoint, lastPoint) - halfWirePitch;
    parameters.m_lineEndZ = std::max(firstPoint, lastPoint) + halfWirePitch;

    return parameters;
  }

} // namespace lar_pandora
