/**
 *  @file   larpandora/LArPandoraInterface/Detectors/LArPandoraDetectorType.h
 *
 *  @brief  Helper functions for extracting detector geometry for use in reconsruction
 *
 *  $Log: $
 */

#ifndef LAR_PANDORA_DETECTOR_TYPE_H
#define LAR_PANDORA_DETECTOR_TYPE_H 1

#include "larcore/Geometry/Geometry.h"

#include "Api/PandoraApi.h"

namespace lar_pandora {

  class LArDriftVolume;
  class LArDetectorGap;
  typedef std::vector<LArDetectorGap> LArDetectorGapList;

  /**
     *  @brief  Empty interface to map pandora to specifics in the LArSoft geometry
     */
  class LArPandoraDetectorType {
  public:
    /**
             *  @brief  Map a LArSoft view to Pandora's U view
             *
             *  @param  tpc the LArSoft TPC ID
             *  @param  cstat the LArSoft cryostat ID
             *  @result The mapped LArSoft view
             */
    virtual geo::View_t TargetViewU(const geo::TPCID::TPCID_t tpc,
                                    const geo::CryostatID::CryostatID_t cstat) const = 0;

    /**
             *  @brief  Map a LArSoft view to Pandora's V view
             *
             *  @param  tpc the LArSoft TPC ID
             *  @param  cstat the LArSoft cryostat ID
             *  @result The mapped LArSoft view
             */
    virtual geo::View_t TargetViewV(const geo::TPCID::TPCID_t tpc,
                                    const geo::CryostatID::CryostatID_t cstat) const = 0;

    /**
             *  @brief  Map a LArSoft view to Pandora's W view
             *
             *  @param  tpc the LArSoft TPC ID
             *  @param  cstat the LArSoft cryostat ID
             *  @result The mapped LArSoft view
             */
    virtual geo::View_t TargetViewW(const geo::TPCID::TPCID_t tpc,
                                    const geo::CryostatID::CryostatID_t cstat) const = 0;

    /**
             *  @brief  The wire pitch of the mapped U view
             *
             *  @result The mapped wire pitch
             */
    virtual float WirePitchU() const = 0;

    /**
             *  @brief  The wire pitch of the mapped V view
             *
             *  @result The mapped wire pitch
             */
    virtual float WirePitchV() const = 0;

    /**
             *  @brief  The wire pitch of the mapped W view
             *
             *  @result The mapped wire pitch
             */
    virtual float WirePitchW() const = 0;

    /**
             *  @brief  The angle of the wires in the mapped U view
             *
             *  @result The mapped wire angle
             */
    virtual float WireAngleU(const geo::TPCID::TPCID_t tpc,
                             const geo::CryostatID::CryostatID_t cstat) const = 0;

    /**
             *  @brief  The angle of the wires in the mapped V view
             *
             *  @result The mapped wire angle
             */
    virtual float WireAngleV(const geo::TPCID::TPCID_t tpc,
                             const geo::CryostatID::CryostatID_t cstat) const = 0;

    /**
             *  @brief  The angle of the wires in the mapped V view
             *
             *  @result The mapped wire angle
             */
    virtual float WireAngleW(const geo::TPCID::TPCID_t tpc,
                             const geo::CryostatID::CryostatID_t cstat) const = 0;

    /**
             *  @brief  Check whether a gap size is small enough to be registered as a detector gap
             *
             *  @param  gaps a cartesean vector holding gap sizes between adjacent TPCs
             *  @param  deltas a cartesean vector holding distances between adjacent TPCs
             *  @param  maxDisplacement the gap size threshold
             *  @result logic bool
             */
    virtual bool CheckDetectorGapSize(const geo::Vector_t& gaps,
                                      const geo::Vector_t& deltas,
                                      const float maxDisplacement) const = 0;

    /**
             *  @brief  Create a detector gap
             *
             *  @param  point1 a point on TPC1 that are closest to TPC2
             *  @param  point2 a point on TPC2 that are closest to TPC1
             *  @param  width the gap sizes
             *  @result the detector gap object
             */
    virtual LArDetectorGap CreateDetectorGap(const geo::Point_t& point1,
                                             const geo::Point_t& point2,
                                             const geo::Vector_t& widths) const = 0;

    /**
             *  @brief  Create detector gaps for all daughter volumes in a logical TPC volume
             *
             *  @param  driftVolume the parent drift volume
             *  @param  maxDisplacement the gap size threshold
             *  @param  listOfGaps the gaps vector to be filled
             */
    virtual void LoadDaughterDetectorGaps(const LArDriftVolume& driftVolume,
                                          const float maxDisplacement,
                                          LArDetectorGapList& listOfGaps) const = 0;

    /**
             *  @brief  Create the line gap parameters to give to the pandora API
             *
             *  @param  gap the input detector gap
             *  @return the pandora API's line gap parameters object
             */
    virtual PandoraApi::Geometry::LineGap::Parameters CreateLineGapParametersFromDetectorGaps(
      const LArDetectorGap& gap) const = 0;

    /**
             *  @brief  Create the line gap parameters to give to the pandora API
             *
             *  @param  view the LArSoft view
             *  @param  tpc the LArSoft TPC ID
             *  @param  cstat the LArSoft cryostat ID
             *  @param  firstXYZ the first 3D coordinate
             *  @param  lastXYZ the last 3D coordinate
             *  @param  halfWirePitch the half wire pitch
             *  @param  pPandora the pandora instance
             *  @param  xFirst the min X of the gap
             *  @param  xLast the max X of the gap
             *  @return the pandora API's line gap parameters object
             */
    virtual PandoraApi::Geometry::LineGap::Parameters CreateLineGapParametersFromReadoutGaps(
      const geo::View_t view,
      const geo::TPCID::TPCID_t tpc,
      const geo::CryostatID::CryostatID_t cstat,
      const double firstXYZ[3],
      const double lastXYZ[3],
      const float halfWirePitch,
      const float xFirst,
      const float xLast,
      const pandora::Pandora* pPandora) const = 0;
  };

  namespace detector_functions {

    /**
         *  @brief  Factory class that returns the correct detector type interface
         *
         *  @result The detector type interface
         */
    LArPandoraDetectorType* GetDetectorType();

    /**
         *  @brief  Calculate the wire angle of a LArTPC view in a given TPC/cryostat
         *
         *  @param  view the LArSoft view
         *  @param  tpc the LArSoft TPC ID
         *  @param  cstat the LArSoft cryostat ID
         *  @param  larsoftGeometry the LArSoft geometry service handle
         *  @return the wire angle
         */
    float WireAngle(const geo::View_t view,
                    const geo::TPCID::TPCID_t tpc,
                    const geo::CryostatID::CryostatID_t cstat,
                    const art::ServiceHandle<geo::Geometry>& larsoftGeometry);

    /**
         *  @brief  Make the drift gap parameters for the Pandora API
         *
         *  @param  gap the detector gap object
         *  @return the line gap parameters for the Pandora API
         */
    PandoraApi::Geometry::LineGap::Parameters CreateDriftGapParameters(const LArDetectorGap& gap);

    /**
         *  @brief  Make the readout gap parameters for the Pandora API
         *
         *  @param  firstPoint the first point in the gap
         *  @param  lastPoint  the last point ih the gap
         *  @param  xFirst the min X of the gap
         *  @param  xLast the max X of the gap
         *  @param  halfWirePitch the half wire pitch
         *  @param  gapType the pandora gap type
         *  @return the line gap parameters for the Pandora API
         */
    PandoraApi::Geometry::LineGap::Parameters CreateReadoutGapParameters(
      const float firstPoint,
      const float lastPoint,
      const float xFirst,
      const float xLast,
      const float halfWirePitch,
      const pandora::LineGapType gapType);

  } // namespace detector_functions

} // namespace lar_pandora
#endif // #ifndef LAR_PANDORA_DETECTOR_TYPE_H
