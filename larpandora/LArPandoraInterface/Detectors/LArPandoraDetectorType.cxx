/**
 *  @file   larpandora/LArPandoraInterface/Detectors/LArPandoraDetectorType.cxx
 *
 *  @brief  Implementation of the interface for handling detector-specific details, as well as some helper functions
 *
 *  $Log: $
 */

#include "larpandora/LArPandoraInterface/Detectors/LArPandoraDetectorType.h"
#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

#include "larpandora/LArPandoraInterface/Detectors/DUNEFarDetVDThreeView.h"
#include "larpandora/LArPandoraInterface/Detectors/ICARUS.h"
#include "larpandora/LArPandoraInterface/Detectors/ProtoDUNEDualPhase.h"
#include "larpandora/LArPandoraInterface/Detectors/VintageLArTPCThreeView.h"

#include <limits>

namespace lar_pandora {

  LArPandoraDetectorType*
  detector_functions::GetDetectorType()
  {
    art::ServiceHandle<geo::Geometry const> geo;

    const unsigned int nPlanes(geo->MaxPlanes());
    std::set<geo::_plane_proj> planeSet;
    for (unsigned int iPlane = 0; iPlane < nPlanes; ++iPlane)
      (void)planeSet.insert(geo->TPC(0, 0).Plane(iPlane).View());

    if (nPlanes == 3 && planeSet.count(geo::kU) && planeSet.count(geo::kY) &&
        planeSet.count(geo::kZ)) {
      return new DUNEFarDetVDThreeView; //TODO Address bare pointer
    }
    else if (nPlanes == 3 && planeSet.count(geo::kU) && planeSet.count(geo::kV) &&
             planeSet.count(geo::kW)) {
      return new VintageLArTPCThreeView;
    }
    else if (nPlanes == 3 && planeSet.count(geo::kU) && planeSet.count(geo::kV) &&
             planeSet.count(geo::kY)) {
      return new ICARUS;
    }
    else if (nPlanes == 2 && planeSet.count(geo::kW) && planeSet.count(geo::kY)) {
      return new ProtoDUNEDualPhase;
    }

    throw cet::exception("LArPandora") << "LArPandoraDetectorType::GetDetectorType --- unable to "
                                          "determine the detector type from the geometry GDML";
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  float
  detector_functions::WireAngle(const geo::View_t view,
                                const geo::TPCID::TPCID_t tpc,
                                const geo::CryostatID::CryostatID_t cstat,
                                const art::ServiceHandle<geo::Geometry>& larsoftGeometry)
  {
    return (0.5f * M_PI - larsoftGeometry->WireAngleToVertical(view, tpc, cstat));
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  PandoraApi::Geometry::LineGap::Parameters
  detector_functions::CreateDriftGapParameters(const LArDetectorGap& gap)
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

  PandoraApi::Geometry::LineGap::Parameters
  detector_functions::CreateReadoutGapParameters(const float firstPoint,
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
