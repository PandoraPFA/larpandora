/**
 *  @file   larpandora/LArPandoraInterface/LArPandoraGeometry.h
 *
 *  @brief  Helper functions for extracting detector geometry for use in reconsruction
 */

#ifndef LAR_PANDORA_GEOMETRY_H
#define LAR_PANDORA_GEOMETRY_H 1

#include "larpandora/LArPandoraInterface/LArPandoraGeometryComponents.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

#include <map>
#include <vector>

namespace lar_pandora {
  class LArPandoraDetectorType;

  /**
 *  @brief  LArPandoraGeometry class
 */
  class LArPandoraGeometry {
  public:
    /**
     *  @brief Load the 2D gaps that go with the chosen geometry
     *
     *  @param listOfGaps the output list of 2D gaps.
     *  @param useActiveBoundingBox when true use ActiveBoundingBox instead of the default midpoint. Meant to handle offsets and things in a better way.
     */
    static void LoadDetectorGaps(LArDetectorGapList& listOfGaps, const bool useActiveBoundingBox);

    /**
     *  @brief Load drift volume geometry
     *
     *  @param outputVolumeList the output list of drift volumes
     *  @param outputVolumeMap the output mapping between cryostat/tpc and drift volumes
     *  @param useActiveBoundingBox when true use ActiveBoundingBox instead of the default midpoint. Meant to handle offsets and things in a better way.
     */
    static void LoadGeometry(LArDriftVolumeList& outputVolumeList,
                             LArDriftVolumeMap& outputVolumeMap,
                             const bool useActiveBoundingBox);

    /**
     *  @brief  Get drift volume ID from a specified cryostat/tpc pair
     *
     *  @param  driftVolumeMap the output mapping between cryostat/tpc and drift volumes
     *  @param  cstat the input cryostat unique ID
     *  @param  tpc the input tpc unique ID
     */
    static unsigned int GetVolumeID(const LArDriftVolumeMap& driftVolumeMap,
                                    const unsigned int cstat,
                                    const unsigned int tpc);

    /**
     *  @brief  Get daughter volume ID from a specified cryostat/tpc pair
     *
     *  @param  driftVolumeMap the output mapping between cryostat/tpc and drift volumes
     *  @param  cstat the input cryostat unique ID
     *  @param  tpc the input tpc unique ID
     */
    static unsigned int GetDaughterVolumeID(const LArDriftVolumeMap& driftVolumeMap,
                                            const unsigned int cstat,
                                            const unsigned int tpc);

    /**
     *  @brief  Convert to global coordinate system
     *
     *  @param  cstat the input cryostat
     *  @param  tpc the input tpc
     *  @param  hit_View the input view
     */
    static geo::View_t GetGlobalView(const unsigned int cstat,
                                     const unsigned int tpc,
                                     const geo::View_t hit_View);

    /**
     *  @brief  Convert to Pandora view
     *
     *  @param  view the input LArSoft view
     *  @param  tpc the input tpc
     *  @param  cstat the input cryostat
     *  @param  detType the input detector type
     *
     *  @return the corresponding Pandora hit type
     */
    static pandora::HitType GetGlobalHitType(const geo::View_t view,
                                             const geo::TPCID::TPCID_t tpc,
                                             const geo::CryostatID::CryostatID_t cstat,
                                             const LArPandoraDetectorType *const detType);

  private:
    /**
     *  @brief  Generate a unique identifier for each TPC
     *
     *  @param  cstat the input cryostat
     *  @param  tpc the input tpc
     */
    static unsigned int GetTpcID(const unsigned int cstat, const unsigned int tpc);

    /**
     *  @brief  Return whether U/V should be switched in global coordinate system for this cryostat/tpc
     *
     *  @param  cstat the input cryostat
     *  @param  tpc the input tpc
     */
    static bool ShouldSwitchUV(const unsigned int cstat, const unsigned int tpc);

    /**
     *  @brief  Return whether U/V should be switched in global coordinate system for this drift direction
     *
     *  @param  isPositiveDrift the drift direction
     */
    static bool ShouldSwitchUV(const bool isPositiveDrift);

    /**
     *  @brief  This method will group TPCs into drift volumes (these are regions of the detector that share a common drift direction,
     *          common range of X coordinates, and common detector parameters such as wire pitch and wire angle).
     *
     *  @param  driftVolumeList to receive the populated drift volume list
     *  @param  useActiveBoundingBox when true use ActiveBoundingBox instead of the default midpoint. Meant to handle offsets and things in a better way.
     */
    static void LoadGeometry(LArDriftVolumeList& driftVolumeList, const bool useActiveBoundingBox);

    /**
     *  @brief  This method will create one or more daughter volumes (these share a common drift orientation along the X-axis,
     *          have parallel or near-parallel wire angles, and similar wire pitches)
     *
     *  @param  driftVolumeList to receive the input drift volume list
     *  @param  parentVolumeList to receive the output daughter drift volume list
     */
    static void LoadGlobalDaughterGeometry(const LArDriftVolumeList& driftVolumeList,
                                           LArDriftVolumeList& daughterVolumeList);

    /**
     *  @brief  This method will return the wire angle for a given hit type, TPC and cryostat
     *
     *  @param  hitType the input hit type
     *  @param  tpc the input TPC
     *  @param  cstat the input cryostat
     *  @param  detType the input detector type
     */
    static float GetWireAngleForHitType(const pandora::HitType hitType, const geo::TPCID::TPCID_t tpc, const geo::CryostatID::CryostatID_t cstat,
                                        const LArPandoraDetectorType *const detType);

    /**
     *  @brief  This method will create one or more readout units, which represent a coherent set of readout planes and their channels (e.g. an APA).
     *
     *  @param  tpcID the input TPC identifier
     *  @param  detType the input detector type
     */
    static LArPandoraReadoutUnitList BuildReadoutUnits(const geo::TPCID &tpcID, const LArPandoraDetectorType *const detType);
  };

} // namespace lar_pandora

#endif // #ifndef LAR_PANDORA_GEOMETRY_H
