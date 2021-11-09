/**
 *  @file   larpandora/LArPandoraInterface/LArPandoraGeometry.cxx
 *
 *  @brief  Helper functions for extracting detector geometry for use in reconsruction
 */

#include "cetlib_except/exception.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

#include "larpandora/LArPandoraInterface/Detectors/LArPandoraDetectorType.h"
#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

#include <iomanip>
#include <set>

namespace lar_pandora {

  void
  LArPandoraGeometry::LoadDetectorGaps(LArDetectorGapList& listOfGaps,
                                       const bool useActiveBoundingBox)
  {
    // Detector gaps can only be loaded once - throw an exception if the output lists are already filled
    if (!listOfGaps.empty())
      throw cet::exception("LArPandora")
        << " LArPandoraGeometry::LoadDetectorGaps --- the list of gaps already exists ";

    // Loop over drift volumes and write out the dead regions at their boundaries
    LArDriftVolumeList driftVolumeList;
    LArPandoraGeometry::LoadGeometry(driftVolumeList, useActiveBoundingBox);

    LArPandoraDetectorType* detType(detector_functions::GetDetectorType());

    for (LArDriftVolumeList::const_iterator iter1 = driftVolumeList.begin(),
                                            iterEnd1 = driftVolumeList.end();
         iter1 != iterEnd1;
         ++iter1) {
      const LArDriftVolume& driftVolume1 = *iter1;

      for (LArDriftVolumeList::const_iterator iter2 = iter1, iterEnd2 = driftVolumeList.end();
           iter2 != iterEnd2;
           ++iter2) {
        const LArDriftVolume& driftVolume2 = *iter2;

        if (driftVolume1.GetVolumeID() == driftVolume2.GetVolumeID()) continue;

        const float maxDisplacement(LArDetectorGap::GetMaxGapSize());

        const float deltaX(std::fabs(driftVolume1.GetCenterX() - driftVolume2.GetCenterX()));
        const float deltaY(std::fabs(driftVolume1.GetCenterY() - driftVolume2.GetCenterY()));
        const float deltaZ(std::fabs(driftVolume1.GetCenterZ() - driftVolume2.GetCenterZ()));

        const float widthX(0.5f * (driftVolume1.GetWidthX() + driftVolume2.GetWidthX()));
        const float widthY(0.5f * (driftVolume1.GetWidthY() + driftVolume2.GetWidthY()));
        const float widthZ(0.5f * (driftVolume1.GetWidthZ() + driftVolume2.GetWidthZ()));

        const float gapX(deltaX - widthX);
        const float gapY(deltaY - widthY);
        const float gapZ(deltaZ - widthZ);

        const float X1((driftVolume1.GetCenterX() < driftVolume2.GetCenterX()) ?
                         (driftVolume1.GetCenterX() + 0.5f * driftVolume1.GetWidthX()) :
                         (driftVolume2.GetCenterX() + 0.5f * driftVolume2.GetWidthX()));
        const float X2((driftVolume1.GetCenterX() > driftVolume2.GetCenterX()) ?
                         (driftVolume1.GetCenterX() - 0.5f * driftVolume1.GetWidthX()) :
                         (driftVolume2.GetCenterX() - 0.5f * driftVolume2.GetWidthX()));
        const float Y1(std::min((driftVolume1.GetCenterY() - 0.5f * driftVolume1.GetWidthY()),
                                (driftVolume2.GetCenterY() - 0.5f * driftVolume2.GetWidthY())));
        const float Y2(std::max((driftVolume1.GetCenterY() + 0.5f * driftVolume1.GetWidthY()),
                                (driftVolume2.GetCenterY() + 0.5f * driftVolume2.GetWidthY())));
        const float Z1(std::min((driftVolume1.GetCenterZ() - 0.5f * driftVolume1.GetWidthZ()),
                                (driftVolume2.GetCenterZ() - 0.5f * driftVolume2.GetWidthZ())));
        const float Z2(std::max((driftVolume1.GetCenterZ() + 0.5f * driftVolume1.GetWidthZ()),
                                (driftVolume2.GetCenterZ() + 0.5f * driftVolume2.GetWidthZ())));

        geo::Vector_t gaps(gapX, gapY, gapZ), deltas(deltaX, deltaY, deltaZ);
        if (detType->CheckDetectorGapSize(gaps, deltas, maxDisplacement)) {
          geo::Point_t point1(X1, Y1, Z1), point2(X2, Y2, Z2);
          geo::Vector_t widths(widthX, widthY, widthZ);
          listOfGaps.emplace_back(detType->CreateDetectorGap(point1, point2, widths));
        }
      }

      detType->LoadDaughterDetectorGaps(driftVolume1, LArDetectorGap::GetMaxGapSize(), listOfGaps);
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void
  LArPandoraGeometry::LoadGeometry(LArDriftVolumeList& outputVolumeList,
                                   LArDriftVolumeMap& outputVolumeMap,
                                   const bool useActiveBoundingBox)
  {
    if (!outputVolumeList.empty())
      throw cet::exception("LArPandora")
        << " LArPandoraGeometry::LoadGeometry --- the list of drift volumes already exists ";

    LArPandoraGeometry::LoadGeometry(outputVolumeList, useActiveBoundingBox);

    // Create mapping between tpc/cstat labels and drift volumes
    for (const LArDriftVolume& driftVolume : outputVolumeList) {
      for (const LArDaughterDriftVolume& tpcVolume : driftVolume.GetTpcVolumeList()) {
        (void)outputVolumeMap.insert(LArDriftVolumeMap::value_type(
          LArPandoraGeometry::GetTpcID(tpcVolume.GetCryostat(), tpcVolume.GetTpc()), driftVolume));
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  unsigned int
  LArPandoraGeometry::GetVolumeID(const LArDriftVolumeMap& driftVolumeMap,
                                  const unsigned int cstat,
                                  const unsigned int tpc)
  {
    if (driftVolumeMap.empty())
      throw cet::exception("LArPandora")
        << " LArPandoraGeometry::GetVolumeID --- detector geometry map is empty";

    LArDriftVolumeMap::const_iterator iter =
      driftVolumeMap.find(LArPandoraGeometry::GetTpcID(cstat, tpc));

    if (driftVolumeMap.end() == iter)
      throw cet::exception("LArPandora")
        << " LArPandoraGeometry::GetVolumeID --- found a TPC that doesn't belong to a drift volume";

    return iter->second.GetVolumeID();
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  unsigned int
  LArPandoraGeometry::GetDaughterVolumeID(const LArDriftVolumeMap& driftVolumeMap,
                                          const unsigned int cstat,
                                          const unsigned int tpc)
  {
    if (driftVolumeMap.empty())
      throw cet::exception("LArPandora")
        << " LArPandoraGeometry::GetDaughterVolumeID --- detector geometry map is empty";

    LArDriftVolumeMap::const_iterator iter =
      driftVolumeMap.find(LArPandoraGeometry::GetTpcID(cstat, tpc));

    if (driftVolumeMap.end() == iter)
      throw cet::exception("LArPandora") << " LArPandoraGeometry::GetDaughterVolumeID --- found a "
                                            "TPC volume that doesn't belong to a drift volume";

    for (LArDaughterDriftVolumeList::const_iterator
           iterDghtr = iter->second.GetTpcVolumeList().begin(),
           iterDghtrEnd = iter->second.GetTpcVolumeList().end();
         iterDghtr != iterDghtrEnd;
         ++iterDghtr) {
      const LArDaughterDriftVolume& daughterVolume(*iterDghtr);
      if (cstat == daughterVolume.GetCryostat() && tpc == daughterVolume.GetTpc())
        return std::distance(iter->second.GetTpcVolumeList().begin(), iterDghtr);
    }
    throw cet::exception("LArPandora")
      << " LArPandoraGeometry::GetDaughterVolumeID --- found a daughter volume that doesn't belong "
         "to the drift volume ";
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  geo::View_t
  LArPandoraGeometry::GetGlobalView(const unsigned int cstat,
                                    const unsigned int tpc,
                                    const geo::View_t hit_View)
  {
    const bool switchUV(LArPandoraGeometry::ShouldSwitchUV(cstat, tpc));

    // ATTN This implicitly assumes that there will be u, v and (maybe) one of either w or y views
    if ((hit_View == geo::kW) || (hit_View == geo::kY)) { return hit_View; }
    else if (hit_View == geo::kU) {
      return (switchUV ? geo::kV : geo::kU);
    }
    else if (hit_View == geo::kV) {
      return (switchUV ? geo::kU : geo::kV);
    }
    else {
      throw cet::exception("LArPandora")
        << " LArPandoraGeometry::GetGlobalView --- found an unknown plane view (not U, V or W) ";
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  unsigned int
  LArPandoraGeometry::GetTpcID(const unsigned int cstat, const unsigned int tpc)
  {
    // We assume there will never be more than 10000 TPCs in a cryostat!
    if (tpc >= 10000)
      throw cet::exception("LArPandora")
        << " LArPandoraGeometry::GetTpcID --- found a TPC with an ID greater than 10000 ";

    return ((10000 * cstat) + tpc);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  bool
  LArPandoraGeometry::ShouldSwitchUV(const unsigned int cstat, const unsigned int tpc)
  {
    // We determine whether U and V views should be switched by checking the drift direction
    art::ServiceHandle<geo::Geometry const> theGeometry;
    const geo::TPCGeo& theTpc(theGeometry->TPC(tpc, cstat));

    const bool isPositiveDrift(theTpc.DriftDirection() == geo::kPosX);
    return LArPandoraGeometry::ShouldSwitchUV(isPositiveDrift);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  bool
  LArPandoraGeometry::ShouldSwitchUV(const bool isPositiveDrift)
  {
    // ATTN: In the dual phase scenario the wire planes pointing along two orthogonal directions and so interchanging U and V is unnecessary
    art::ServiceHandle<geo::Geometry const> theGeometry;
    if (theGeometry->MaxPlanes() == 2) return false;

    // We assume that all multiple drift volume detectors have the APA - CPA - APA - CPA design
    return isPositiveDrift;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void
  LArPandoraGeometry::LoadGeometry(LArDriftVolumeList& driftVolumeList,
                                   const bool useActiveBoundingBox)
  {
    // This method will group TPCs into "drift volumes" (these are regions of the detector that share a common drift direction,
    // common range of x coordinates, and common detector parameters such as wire pitch and wire angle).
    if (!driftVolumeList.empty())
      throw cet::exception("LArPandora")
        << " LArPandoraGeometry::LoadGeometry --- detector geometry has already been loaded ";

    typedef std::set<unsigned int> UIntSet;

    // Pandora requires three independent images, and ability to correlate features between images (via wire angles and transformation plugin).
    art::ServiceHandle<geo::Geometry const> theGeometry;
    LArPandoraDetectorType* detType(detector_functions::GetDetectorType());
    const float wirePitchU(detType->WirePitchU());
    const float wirePitchV(detType->WirePitchV());
    const float wirePitchW(detType->WirePitchW());
    const float maxDeltaTheta(0.01f); // leave this hard-coded for now

    // Loop over cryostats
    for (unsigned int icstat = 0; icstat < theGeometry->Ncryostats(); ++icstat) {
      UIntSet cstatList;

      // Loop over TPCs in in this cryostat
      for (unsigned int itpc1 = 0; itpc1 < theGeometry->NTPC(icstat); ++itpc1) {
        if (cstatList.end() != cstatList.find(itpc1)) continue;

        // Use this TPC to seed a drift volume
        const geo::TPCGeo& theTpc1(theGeometry->TPC(itpc1, icstat));
        cstatList.insert(itpc1);

        const float wireAngleU(detType->WireAngleU(itpc1, icstat));
        const float wireAngleV(detType->WireAngleV(itpc1, icstat));
        const float wireAngleW(detType->WireAngleW(itpc1, icstat));

        double localCoord1[3] = {0., 0., 0.};
        double worldCoord1[3] = {0., 0., 0.};
        theTpc1.LocalToWorld(localCoord1, worldCoord1);

        float driftMinX(useActiveBoundingBox ? theTpc1.ActiveBoundingBox().MinX() :
                                               (worldCoord1[0] - theTpc1.ActiveHalfWidth()));
        float driftMaxX(useActiveBoundingBox ? theTpc1.ActiveBoundingBox().MaxX() :
                                               (worldCoord1[0] + theTpc1.ActiveHalfWidth()));
        float driftMinY(useActiveBoundingBox ? theTpc1.ActiveBoundingBox().MinY() :
                                               (worldCoord1[1] - theTpc1.ActiveHalfHeight()));
        float driftMaxY(useActiveBoundingBox ? theTpc1.ActiveBoundingBox().MaxY() :
                                               (worldCoord1[1] + theTpc1.ActiveHalfHeight()));
        float driftMinZ(useActiveBoundingBox ? theTpc1.ActiveBoundingBox().MinZ() :
                                               (worldCoord1[2] - 0.5f * theTpc1.ActiveLength()));
        float driftMaxZ(useActiveBoundingBox ? theTpc1.ActiveBoundingBox().MaxZ() :
                                               (worldCoord1[2] + 0.5f * theTpc1.ActiveLength()));

        const double min1(
          useActiveBoundingBox ?
            (0.5 * (driftMinX + driftMaxX) - 0.25 * std::fabs(driftMaxX - driftMinX)) :
            (worldCoord1[0] - 0.5 * theTpc1.ActiveHalfWidth()));
        const double max1(
          useActiveBoundingBox ?
            (0.5 * (driftMinX + driftMaxX) + 0.25 * std::fabs(driftMaxX - driftMinX)) :
            (worldCoord1[0] + 0.5 * theTpc1.ActiveHalfWidth()));

        const bool isPositiveDrift(theTpc1.DriftDirection() == geo::kPosX);

        UIntSet tpcList;
        tpcList.insert(itpc1);

        LArDaughterDriftVolumeList tpcVolumeList;
        tpcVolumeList.emplace_back(LArDaughterDriftVolume(icstat,
                                                          itpc1,
                                                          0.5f * (driftMaxX + driftMinX),
                                                          0.5f * (driftMaxY + driftMinY),
                                                          0.5f * (driftMaxZ + driftMinZ),
                                                          (driftMaxX - driftMinX),
                                                          (driftMaxY - driftMinY),
                                                          (driftMaxZ - driftMinZ)));

        // Now identify the other TPCs associated with this drift volume
        for (unsigned int itpc2 = itpc1 + 1; itpc2 < theGeometry->NTPC(icstat); ++itpc2) {
          if (cstatList.end() != cstatList.find(itpc2)) continue;

          const geo::TPCGeo& theTpc2(theGeometry->TPC(itpc2, icstat));

          if (theTpc1.DriftDirection() != theTpc2.DriftDirection()) continue;

          const float dThetaU(detType->WireAngleU(itpc1, icstat) -
                              detType->WireAngleU(itpc2, icstat));
          const float dThetaV(detType->WireAngleV(itpc1, icstat) -
                              detType->WireAngleV(itpc2, icstat));
          const float dThetaW(detType->WireAngleW(itpc1, icstat) -
                              detType->WireAngleW(itpc2, icstat));
          if (dThetaU > maxDeltaTheta || dThetaV > maxDeltaTheta || dThetaW > maxDeltaTheta)
            continue;

          double localCoord2[3] = {0., 0., 0.};
          double worldCoord2[3] = {0., 0., 0.};
          theTpc2.LocalToWorld(localCoord2, worldCoord2);

          const float driftMinX2(useActiveBoundingBox ?
                                   theTpc2.ActiveBoundingBox().MinX() :
                                   (worldCoord2[0] - theTpc2.ActiveHalfWidth()));
          const float driftMaxX2(useActiveBoundingBox ?
                                   theTpc2.ActiveBoundingBox().MaxX() :
                                   (worldCoord2[0] + theTpc2.ActiveHalfWidth()));

          const double min2(
            useActiveBoundingBox ?
              (0.5 * (driftMinX2 + driftMaxX2) - 0.25 * std::fabs(driftMaxX2 - driftMinX2)) :
              (worldCoord2[0] - 0.5 * theTpc2.ActiveHalfWidth()));
          const double max2(
            useActiveBoundingBox ?
              (0.5 * (driftMinX2 + driftMaxX2) + 0.25 * std::fabs(driftMaxX2 - driftMinX2)) :
              (worldCoord2[0] + 0.5 * theTpc2.ActiveHalfWidth()));

          if ((min2 > max1) || (min1 > max2)) continue;

          cstatList.insert(itpc2);
          tpcList.insert(itpc2);

          const float driftMinY2(useActiveBoundingBox ?
                                   theTpc2.ActiveBoundingBox().MinY() :
                                   (worldCoord2[1] - theTpc2.ActiveHalfHeight()));
          const float driftMaxY2(useActiveBoundingBox ?
                                   theTpc2.ActiveBoundingBox().MaxY() :
                                   (worldCoord2[1] + theTpc2.ActiveHalfHeight()));
          const float driftMinZ2(useActiveBoundingBox ?
                                   theTpc2.ActiveBoundingBox().MinZ() :
                                   (worldCoord2[2] - 0.5f * theTpc2.ActiveLength()));
          const float driftMaxZ2(useActiveBoundingBox ?
                                   theTpc2.ActiveBoundingBox().MaxZ() :
                                   (worldCoord2[2] + 0.5f * theTpc2.ActiveLength()));

          driftMinX = std::min(driftMinX, driftMinX2);
          driftMaxX = std::max(driftMaxX, driftMaxX2);
          driftMinY = std::min(driftMinY, driftMinY2);
          driftMaxY = std::max(driftMaxY, driftMaxY2);
          driftMinZ = std::min(driftMinZ, driftMinZ2);
          driftMaxZ = std::max(driftMaxZ, driftMaxZ2);

          tpcVolumeList.emplace_back(LArDaughterDriftVolume(icstat,
                                                            itpc2,
                                                            0.5f * (driftMaxX2 + driftMinX2),
                                                            0.5f * (driftMaxY2 + driftMinY2),
                                                            0.5f * (driftMaxZ2 + driftMinZ2),
                                                            (driftMaxX2 - driftMinX2),
                                                            (driftMaxY2 - driftMinY2),
                                                            (driftMaxZ2 - driftMinZ2)));
        }

        // Create new daughter drift volume (volume ID = 0 to N-1)
        driftVolumeList.emplace_back(driftVolumeList.size(),
                                     isPositiveDrift,
                                     wirePitchU,
                                     wirePitchV,
                                     wirePitchW,
                                     wireAngleU,
                                     wireAngleV,
                                     wireAngleW,
                                     0.5f * (driftMaxX + driftMinX),
                                     0.5f * (driftMaxY + driftMinY),
                                     0.5f * (driftMaxZ + driftMinZ),
                                     (driftMaxX - driftMinX),
                                     (driftMaxY - driftMinY),
                                     (driftMaxZ - driftMinZ),
                                     (wirePitchU + wirePitchV + wirePitchW + 0.1f),
                                     tpcVolumeList);
      }
    }

    if (driftVolumeList.empty())
      throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadGeometry --- failed to find "
                                            "any drift volumes in this detector geometry ";
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void
  LArPandoraGeometry::LoadGlobalDaughterGeometry(const LArDriftVolumeList& driftVolumeList,
                                                 LArDriftVolumeList& daughterVolumeList)
  {
    // This method will create one or more daughter volumes (these share a common drift orientation along the X-axis,
    // have parallel or near-parallel wire angles, and similar wire pitches)
    //
    // ATTN: we assume that the U and V planes have equal and opposite wire orientations

    if (!daughterVolumeList.empty())
      throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadGlobalDaughterGeometry --- "
                                            "daughter geometry has already been loaded ";

    if (driftVolumeList.empty())
      throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadGlobalDaughterGeometry --- "
                                            "detector geometry has not yet been loaded ";

    std::cout << "The size of the drif list is: " << driftVolumeList.size() << std::endl;
    int count(0);
    // Create daughter drift volumes
    for (const LArDriftVolume& driftVolume : driftVolumeList) {
      std::cout << "Looking at dau vol: " << count++ << std::endl;
      const bool switchViews(LArPandoraGeometry::ShouldSwitchUV(driftVolume.IsPositiveDrift()));

      const float daughterWirePitchU(switchViews ? driftVolume.GetWirePitchV() :
                                                   driftVolume.GetWirePitchU());
      const float daughterWirePitchV(switchViews ? driftVolume.GetWirePitchU() :
                                                   driftVolume.GetWirePitchV());
      const float daughterWirePitchW(driftVolume.GetWirePitchW());
      const float daughterWireAngleU(switchViews ? driftVolume.GetWireAngleV() :
                                                   driftVolume.GetWireAngleU());
      const float daughterWireAngleV(switchViews ? driftVolume.GetWireAngleU() :
                                                   driftVolume.GetWireAngleV());
      const float daughterWireAngleW(driftVolume.GetWireAngleW());

      daughterVolumeList.push_back(LArDriftVolume(driftVolume.GetVolumeID(),
                                                  driftVolume.IsPositiveDrift(),
                                                  daughterWirePitchU,
                                                  daughterWirePitchV,
                                                  daughterWirePitchW,
                                                  daughterWireAngleU,
                                                  daughterWireAngleV,
                                                  daughterWireAngleW,
                                                  driftVolume.GetCenterX(),
                                                  driftVolume.GetCenterY(),
                                                  driftVolume.GetCenterZ(),
                                                  driftVolume.GetWidthX(),
                                                  driftVolume.GetWidthY(),
                                                  driftVolume.GetWidthZ(),
                                                  driftVolume.GetSigmaUVZ(),
                                                  driftVolume.GetTpcVolumeList()));
    }

    if (daughterVolumeList.empty())
      throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadGlobalDaughterGeometry --- "
                                            "failed to create daughter geometry list ";
  }

  //------------------------------------------------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------------------------------------------------

  LArDriftVolume::LArDriftVolume(const unsigned int volumeID,
                                 const bool isPositiveDrift,
                                 const float wirePitchU,
                                 const float wirePitchV,
                                 const float wirePitchW,
                                 const float wireAngleU,
                                 const float wireAngleV,
                                 const float wireAngleW,
                                 const float centerX,
                                 const float centerY,
                                 const float centerZ,
                                 const float widthX,
                                 const float widthY,
                                 const float widthZ,
                                 const float sigmaUVZ,
                                 const LArDaughterDriftVolumeList& tpcVolumeList)
    : m_volumeID(volumeID)
    , m_isPositiveDrift(isPositiveDrift)
    , m_wirePitchU(wirePitchU)
    , m_wirePitchV(wirePitchV)
    , m_wirePitchW(wirePitchW)
    , m_wireAngleU(wireAngleU)
    , m_wireAngleV(wireAngleV)
    , m_wireAngleW(wireAngleW)
    , m_centerX(centerX)
    , m_centerY(centerY)
    , m_centerZ(centerZ)
    , m_widthX(widthX)
    , m_widthY(widthY)
    , m_widthZ(widthZ)
    , m_sigmaUVZ(sigmaUVZ)
    , m_tpcVolumeList(tpcVolumeList)
  {}

} // namespace lar_pandora
