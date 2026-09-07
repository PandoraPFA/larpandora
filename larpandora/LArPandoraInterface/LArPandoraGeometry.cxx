/**
 *  @file   larpandora/LArPandoraInterface/LArPandoraGeometry.cxx
 *
 *  @brief  Helper functions for extracting detector geometry for use in reconsruction
 */

#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib_except/exception.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

#include "larpandora/LArPandoraInterface/Detectors/GetDetectorType.h"
#include "larpandora/LArPandoraInterface/Detectors/LArPandoraDetectorType.h"

#include "larpandoracontent/LArPersistency/LArTPCFactory.h"

#include <iomanip>
#include <set>

namespace lar_pandora {

  void LArPandoraGeometry::LoadDetectorGaps(LArDetectorGapList& listOfGaps,
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

  void LArPandoraGeometry::LoadGeometry(LArDriftVolumeList& outputVolumeList,
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

  unsigned int LArPandoraGeometry::GetVolumeID(const LArDriftVolumeMap& driftVolumeMap,
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

  unsigned int LArPandoraGeometry::GetDaughterVolumeID(const LArDriftVolumeMap& driftVolumeMap,
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

  geo::View_t LArPandoraGeometry::GetGlobalView(const unsigned int cstat,
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

    throw cet::exception("LArPandora")
      << " LArPandoraGeometry::GetGlobalView --- found an unknown plane view (not U, V or W) ";
  }

  pandora::HitType LArPandoraGeometry::GetGlobalHitType(const geo::View_t view,
                                                        const geo::TPCID::TPCID_t tpc,
                                                        const geo::CryostatID::CryostatID_t cstat,
                                                        const LArPandoraDetectorType* const detType)
  {
    if (view == detType->TargetViewW(tpc, cstat)) return pandora::TPC_VIEW_W;
    if (view == detType->TargetViewU(tpc, cstat)) return pandora::TPC_VIEW_U;
    if (view == detType->TargetViewV(tpc, cstat)) return pandora::TPC_VIEW_V;
    throw cet::exception("LArPandora") << " GetGlobalHitType --- unrecognised wire view ";
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  unsigned int LArPandoraGeometry::GetTpcID(const unsigned int cstat, const unsigned int tpc)
  {
    // We assume there will never be more than 10000 TPCs in a cryostat!
    if (tpc >= 10000)
      throw cet::exception("LArPandora")
        << " LArPandoraGeometry::GetTpcID --- found a TPC with an ID greater than 10000 ";

    return ((10000 * cstat) + tpc);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  bool LArPandoraGeometry::ShouldSwitchUV(const unsigned int cstat, const unsigned int tpc)
  {
    // We determine whether U and V views should be switched by checking the drift direction
    art::ServiceHandle<geo::Geometry const> theGeometry;

    const bool isPositiveDrift(theGeometry->TPC({cstat, tpc}).DriftSign() ==
                               geo::DriftSign::Positive);
    return LArPandoraGeometry::ShouldSwitchUV(isPositiveDrift);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  bool LArPandoraGeometry::ShouldSwitchUV(const bool isPositiveDrift)
  {
    // ATTN: In the dual phase scenario the wire planes pointing along two orthogonal directions and so interchanging U and V is unnecessary
    auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout>()->Get();
    if (wireReadoutGeom.MaxPlanes() == 2) return false;

    // We assume that all multiple drift volume detectors have the APA - CPA - APA - CPA design
    return isPositiveDrift;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void LArPandoraGeometry::LoadGeometry(LArDriftVolumeList& driftVolumeList,
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
    for (auto const& cryostat : theGeometry->Iterate<geo::CryostatGeo>()) {
      auto const icstat = cryostat.ID().Cryostat;
      UIntSet cstatList;

      // Loop over TPCs in in this cryostat
      for (auto const& theTpc1 : theGeometry->Iterate<geo::TPCGeo>(cryostat.ID())) {
        auto const itpc1 = theTpc1.ID().TPC;
        if (cstatList.end() != cstatList.find(itpc1)) continue;

        // Use this TPC to seed a drift volume
        cstatList.insert(itpc1);

        const float wireAngleU(detType->WireAngleU(itpc1, icstat));
        const float wireAngleV(detType->WireAngleV(itpc1, icstat));
        const float wireAngleW(detType->WireAngleW(itpc1, icstat));

        auto const worldCoord1 = theTpc1.GetCenter();

        float driftMinX(useActiveBoundingBox ? theTpc1.ActiveBoundingBox().MinX() :
                                               (worldCoord1.X() - theTpc1.ActiveHalfWidth()));
        float driftMaxX(useActiveBoundingBox ? theTpc1.ActiveBoundingBox().MaxX() :
                                               (worldCoord1.X() + theTpc1.ActiveHalfWidth()));
        float driftMinY(useActiveBoundingBox ? theTpc1.ActiveBoundingBox().MinY() :
                                               (worldCoord1.Y() - theTpc1.ActiveHalfHeight()));
        float driftMaxY(useActiveBoundingBox ? theTpc1.ActiveBoundingBox().MaxY() :
                                               (worldCoord1.Y() + theTpc1.ActiveHalfHeight()));
        float driftMinZ(useActiveBoundingBox ? theTpc1.ActiveBoundingBox().MinZ() :
                                               (worldCoord1.Z() - 0.5f * theTpc1.ActiveLength()));
        float driftMaxZ(useActiveBoundingBox ? theTpc1.ActiveBoundingBox().MaxZ() :
                                               (worldCoord1.Z() + 0.5f * theTpc1.ActiveLength()));

        const double min1(
          useActiveBoundingBox ?
            (0.5 * (driftMinX + driftMaxX) - 0.25 * std::fabs(driftMaxX - driftMinX)) :
            (worldCoord1.X() - 0.5 * theTpc1.ActiveHalfWidth()));
        const double max1(
          useActiveBoundingBox ?
            (0.5 * (driftMinX + driftMaxX) + 0.25 * std::fabs(driftMaxX - driftMinX)) :
            (worldCoord1.X() + 0.5 * theTpc1.ActiveHalfWidth()));

        const bool isPositiveDrift(theTpc1.DriftSign() == geo::DriftSign::Positive);

        UIntSet tpcList;
        tpcList.insert(itpc1);

        LArDaughterDriftVolumeList tpcVolumeList;
        tpcVolumeList.emplace_back(
          LArDaughterDriftVolume(icstat,
                                 itpc1,
                                 0.5f * (driftMaxX + driftMinX),
                                 0.5f * (driftMaxY + driftMinY),
                                 0.5f * (driftMaxZ + driftMinZ),
                                 (driftMaxX - driftMinX),
                                 (driftMaxY - driftMinY),
                                 (driftMaxZ - driftMinZ),
                                 BuildReadoutUnits(geo::TPCID{icstat, itpc1}, detType)));

        // Now identify the other TPCs associated with this drift volume
        for (auto const& theTpc2 : theGeometry->Iterate<geo::TPCGeo>(cryostat.ID())) {
          auto const itpc2 = theTpc2.ID().TPC;
          if (cstatList.end() != cstatList.find(itpc2)) continue;
          if (theTpc1.DriftSign() != theTpc2.DriftSign()) continue;

          const float dThetaU(detType->WireAngleU(itpc1, icstat) -
                              detType->WireAngleU(itpc2, icstat));
          const float dThetaV(detType->WireAngleV(itpc1, icstat) -
                              detType->WireAngleV(itpc2, icstat));
          const float dThetaW(detType->WireAngleW(itpc1, icstat) -
                              detType->WireAngleW(itpc2, icstat));
          if (dThetaU > maxDeltaTheta || dThetaV > maxDeltaTheta || dThetaW > maxDeltaTheta)
            continue;

          auto const worldCoord2 = theTpc2.GetCenter();

          const float driftMinX2(useActiveBoundingBox ?
                                   theTpc2.ActiveBoundingBox().MinX() :
                                   (worldCoord2.X() - theTpc2.ActiveHalfWidth()));
          const float driftMaxX2(useActiveBoundingBox ?
                                   theTpc2.ActiveBoundingBox().MaxX() :
                                   (worldCoord2.X() + theTpc2.ActiveHalfWidth()));

          const double min2(
            useActiveBoundingBox ?
              (0.5 * (driftMinX2 + driftMaxX2) - 0.25 * std::fabs(driftMaxX2 - driftMinX2)) :
              (worldCoord2.X() - 0.5 * theTpc2.ActiveHalfWidth()));
          const double max2(
            useActiveBoundingBox ?
              (0.5 * (driftMinX2 + driftMaxX2) + 0.25 * std::fabs(driftMaxX2 - driftMinX2)) :
              (worldCoord2.X() + 0.5 * theTpc2.ActiveHalfWidth()));

          if ((min2 > max1) || (min1 > max2)) continue;

          cstatList.insert(itpc2);
          tpcList.insert(itpc2);

          const float driftMinY2(useActiveBoundingBox ?
                                   theTpc2.ActiveBoundingBox().MinY() :
                                   (worldCoord2.Y() - theTpc2.ActiveHalfHeight()));
          const float driftMaxY2(useActiveBoundingBox ?
                                   theTpc2.ActiveBoundingBox().MaxY() :
                                   (worldCoord2.Y() + theTpc2.ActiveHalfHeight()));
          const float driftMinZ2(useActiveBoundingBox ?
                                   theTpc2.ActiveBoundingBox().MinZ() :
                                   (worldCoord2.Z() - 0.5f * theTpc2.ActiveLength()));
          const float driftMaxZ2(useActiveBoundingBox ?
                                   theTpc2.ActiveBoundingBox().MaxZ() :
                                   (worldCoord2.Z() + 0.5f * theTpc2.ActiveLength()));

          driftMinX = std::min(driftMinX, driftMinX2);
          driftMaxX = std::max(driftMaxX, driftMaxX2);
          driftMinY = std::min(driftMinY, driftMinY2);
          driftMaxY = std::max(driftMaxY, driftMaxY2);
          driftMinZ = std::min(driftMinZ, driftMinZ2);
          driftMaxZ = std::max(driftMaxZ, driftMaxZ2);

          tpcVolumeList.emplace_back(
            LArDaughterDriftVolume(icstat,
                                   itpc2,
                                   0.5f * (driftMaxX2 + driftMinX2),
                                   0.5f * (driftMaxY2 + driftMinY2),
                                   0.5f * (driftMaxZ2 + driftMinZ2),
                                   (driftMaxX2 - driftMinX2),
                                   (driftMaxY2 - driftMinY2),
                                   (driftMaxZ2 - driftMinZ2),
                                   BuildReadoutUnits(geo::TPCID{icstat, itpc2}, detType)));
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

  void LArPandoraGeometry::LoadGlobalDaughterGeometry(const LArDriftVolumeList& driftVolumeList,
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

  float LArPandoraGeometry::GetWireAngleForHitType(const pandora::HitType hitType,
                                                   const geo::TPCID::TPCID_t tpc,
                                                   const geo::CryostatID::CryostatID_t cstat,
                                                   const LArPandoraDetectorType* const detType)
  {
    if (pandora::TPC_VIEW_U == hitType) return detType->WireAngleU(tpc, cstat);
    if (pandora::TPC_VIEW_V == hitType) return detType->WireAngleV(tpc, cstat);
    if (pandora::TPC_VIEW_W == hitType) return detType->WireAngleW(tpc, cstat);

    throw cet::exception("LArPandora") << " GetWireAngleForHitType --- unrecognised hit type ";
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  float LArPandoraGeometry::GetWirePitchForHitType(const pandora::HitType hitType,
                                                   const LArPandoraDetectorType* const detType)
  {
    if (pandora::TPC_VIEW_U == hitType) return detType->WirePitchU();
    if (pandora::TPC_VIEW_V == hitType) return detType->WirePitchV();
    if (pandora::TPC_VIEW_W == hitType) return detType->WirePitchW();

    throw cet::exception("LArPandora") << " GetWirePitchForHitType --- unrecognised hit type ";
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  LArPandoraReadoutUnitList LArPandoraGeometry::BuildReadoutUnits(
    const geo::TPCID& tpcID,
    const LArPandoraDetectorType* const detType)
  {
    auto const& channelReadout{art::ServiceHandle<geo::WireReadout const>()->Get()};
    LArPandoraReadoutUnitList readoutUnitList;

    // First pass: map Pandora views to planes (via detType - handles HD's U/V swap).
    std::map<pandora::HitType, const geo::PlaneGeo*> hitTypeToPlane;
    for (auto const& plane : channelReadout.Iterate<geo::PlaneGeo>(tpcID)) {
      const pandora::HitType globalHitType(
        LArPandoraGeometry::GetGlobalHitType(plane.View(), tpcID.TPC, tpcID.Cryostat, detType));
      hitTypeToPlane[globalHitType] = &plane;
    }

    // Second pass: compute the generative parameters used to compute intervals.
    struct PlaneInfo {
      float m_angle;
      float m_referenceCoordinate;
      float m_pitch;
      pandora::CartesianVector m_unitCenter{0.f, 0.f, 0.f};
      pandora::CartesianVector m_unitSize{0.f, 0.f, 0.f};
      unsigned int m_nChannels;
    };

    std::map<pandora::HitType, PlaneInfo> hitTypeToPlaneInfo;
    for (auto const& [hitType, pPlane] : hitTypeToPlane) {
      PlaneInfo info;
      info.m_angle =
        LArPandoraGeometry::GetWireAngleForHitType(hitType, tpcID.TPC, tpcID.Cryostat, detType);

      const geo::WireGeo& wire0(pPlane->Wire(0));
      const geo::WireGeo& wire1(pPlane->Wire(1));
      const auto projectCoordinate = [&info](const geo::Point_t& p) {
        return p.Z() * std::cos(info.m_angle) - p.Y() * std::sin(info.m_angle);
      };

      info.m_referenceCoordinate = projectCoordinate(wire0.GetCenter());
      const float rawPitch(projectCoordinate(wire1.GetCenter()) - info.m_referenceCoordinate);
      const float nominalPitch(LArPandoraGeometry::GetWirePitchForHitType(hitType, detType));
      info.m_pitch = std::copysign(nominalPitch, rawPitch);

      const geo::BoxBoundedGeo box(pPlane->BoundingBox());
      info.m_unitCenter = pandora::CartesianVector(
        0.f, 0.5f * (box.MinY() + box.MaxY()), 0.5f * (box.MinZ() + box.MaxZ()));
      info.m_unitSize =
        pandora::CartesianVector(0.f, box.MaxY() - box.MinY(), box.MaxZ() - box.MinZ());
      info.m_nChannels = channelReadout.Nwires(pPlane->ID());

      hitTypeToPlaneInfo[hitType] = info;
    }

    // Third pass: compute the channel overlap intervals, via the same ComputeChannelInterval logic used to regenerate this data from a
    // persisted geometry file when running standalone.
    for (auto const& [hitType, pPlane] : hitTypeToPlane) {
      const PlaneInfo& selfInfo(hitTypeToPlaneInfo.at(hitType));

      LArPandoraReadoutChannelList channelList;
      channelList.reserve(selfInfo.m_nChannels);
      for (unsigned int iChannel = 0; iChannel < selfInfo.m_nChannels; ++iChannel) {
        const float selfCoordinate(selfInfo.m_referenceCoordinate +
                                   static_cast<float>(iChannel) * selfInfo.m_pitch);

        pandora::LArReadoutChannel::ViewChannelIntervalArray intervals;
        std::size_t slot(0);

        for (auto const& [otherHitType, otherInfo] : hitTypeToPlaneInfo) {
          if (otherHitType == hitType) continue;

          intervals[slot++] = {
            otherHitType,
            lar_content::LArTPCFactory::ComputeChannelInterval(selfCoordinate,
                                                               selfInfo.m_angle,
                                                               selfInfo.m_unitCenter,
                                                               selfInfo.m_unitSize,
                                                               otherInfo.m_angle,
                                                               otherInfo.m_referenceCoordinate,
                                                               otherInfo.m_pitch,
                                                               otherInfo.m_nChannels)};
        }

        channelList.emplace_back(iChannel, intervals);
      }

      readoutUnitList.emplace_back(pPlane->ID().Plane,
                                   hitType,
                                   selfInfo.m_referenceCoordinate,
                                   selfInfo.m_pitch,
                                   selfInfo.m_unitCenter,
                                   selfInfo.m_unitSize,
                                   channelList);
    }

    return readoutUnitList;
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
