/**
 *  @file   larpandora/LArPandoraInterface/Detectors/DUNEFarDetVDThreeView30DegDriftY.h
 *
 *  @brief  Detector interface for DUNE's vertical drift, 3 view far detector
 *
 *  $Log: $
 */

#include "larpandora/LArPandoraInterface/Detectors/VintageLArTPCThreeView.h"

#include "larcore/Geometry/Geometry.h"

namespace lar_pandora {

  /**
     *  @brief  Detector interface DUNE's vertical drift 30 degree view far detector with Y as the drift
     */
  class DUNEFarDetVDThreeView30DegDriftY : public VintageLArTPCThreeView {
  public:
      DUNEFarDetVDThreeView30DegDriftY();

      geo::Point_t RotateToDriftX(const geo::Point_t& globalPoint) const override;

      geo::Point_t UndoRotateToDriftX(const geo::Point_t& localPoint) const override;

  private:
      geo::Rotation_t m_DriftXRotation;        ///< the rotation to make the drift axis the x-axis
      geo::Rotation_t m_InverseDriftXRotation; ///< the rotation that undoes the drift==x rotation
  };

  //------------------------------------------------------------------------------------------------------------------------------------------

  DUNEFarDetVDThreeView30DegDriftY::DUNEFarDetVDThreeView30DegDriftY() :
      VintageLArTPCThreeView(),
      m_DriftXRotation (0., -1., 0.,
                        1.,  0., 0.,
                        0.,  0., 1.),
      m_InverseDriftXRotation ( 0.,  1., 0.,
                               -1.,  0., 0.,
                                0.,  0., 1.) {}


  //------------------------------------------------------------------------------------------------------------------------------------------

  geo::Point_t DUNEFarDetVDThreeView30DegDriftY::RotateToDriftX(const geo::Point_t& globalPoint) const
  {
      return m_DriftXRotation*globalPoint;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  geo::Point_t DUNEFarDetVDThreeView30DegDriftY::UndoRotateToDriftX(const geo::Point_t& localPoint) const
  {
      return m_InverseDriftXRotation*localPoint;
  }

} // namespace lar_pandora
