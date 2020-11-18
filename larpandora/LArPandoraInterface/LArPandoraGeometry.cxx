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

#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <iomanip>
#include <set>

namespace lar_pandora {

  void
  LArPandoraGeometry::LoadDetectorGaps(LArDetectorGapList& listOfGaps)
  {
    // Detector gaps can only be loaded once - throw an exception if the output lists are already filled
    if (!listOfGaps.empty())
      throw cet::exception("LArPandora")
        << " LArPandoraGeometry::LoadDetectorGaps --- the list of gaps already exists ";

    // Loop over drift volumes and write out the dead regions at their boundaries
    LArDriftVolumeList driftVolumeList;
    LArPandoraGeometry::LoadGeometry(driftVolumeList);

    // ATTN: Expectations here are that the input geometry corresponds to either a single or dual phase LArTPC.
    art::ServiceHandle<geo::Geometry const> theGeometry;
    const bool isDualPhase(LArPandoraHelper::IsDualPhase());

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

        if (!isDualPhase && (gapX < 0.f || gapX > maxDisplacement || deltaY > maxDisplacement ||
                             deltaZ > maxDisplacement))
          continue;

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

        if (isDualPhase && (std::fabs(gapY) > maxDisplacement || std::fabs(gapZ) > maxDisplacement))
          listOfGaps.emplace_back(
            LArDetectorGap(X1, Y1 + widthY, Z1 + widthZ, X2, Y2 - widthY, Z2 - widthZ));

        else if (!isDualPhase)
          listOfGaps.emplace_back(LArDetectorGap(X1, Y1, Z1, X2, Y2, Z2));
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void
  LArPandoraGeometry::LoadGeometry(LArDriftVolumeList& outputVolumeList,
                                   LArDriftVolumeMap& outputVolumeMap)
  {
    if (!outputVolumeList.empty())
      throw cet::exception("LArPandora")
        << " LArPandoraGeometry::LoadGeometry --- the list of drift volumes already exists ";

    // Use a global coordinate system but keep drift volumes separate
    LArDriftVolumeList inputVolumeList;
    LArPandoraGeometry::LoadGeometry(inputVolumeList);
    LArPandoraGeometry::LoadGlobalDaughterGeometry(inputVolumeList, outputVolumeList);

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
  LArPandoraGeometry::LoadGeometry(LArDriftVolumeList& driftVolumeList)
  {
    // This method will group TPCs into "drift volumes" (these are regions of the detector that share a common drift direction,
    // common range of x coordinates, and common detector parameters such as wire pitch and wire angle).
    if (!driftVolumeList.empty())
      throw cet::exception("LArPandora")
        << " LArPandoraGeometry::LoadGeometry --- detector geometry has already been loaded ";

    typedef std::set<unsigned int> UIntSet;

    // Pandora requires three independent images, and ability to correlate features between images (via wire angles and transformation plugin).
    art::ServiceHandle<geo::Geometry const> theGeometry;
    const unsigned int nWirePlanes(theGeometry->MaxPlanes());

    if (nWirePlanes > 3)
      throw cet::exception("LArPandora")
        << " LArPandoraGeometry::LoadGeometry --- More than three wire planes present ";

    // We here check the plane information only for the first tpc in the first cryostat.
    if ((0 == theGeometry->Ncryostats()) || (0 == theGeometry->NTPC(0)))
      throw cet::exception("LArPandora")
        << " LArPandoraGeometry::LoadGeometry --- unable to access first tpc in first cryostat ";

    std::unordered_set<geo::_plane_proj> planeSet;
    for (unsigned int iPlane = 0; iPlane < nWirePlanes; ++iPlane)
      (void)planeSet.insert(theGeometry->TPC(0, 0).Plane(iPlane).View());

    // ATTN: Expectations here are that the input geometry corresponds to either a single or dual phase LArTPC.  For single phase we expect
    // three views, U, V and either W or Y, for dual phase we expect two views, W and Y.
    const bool isDualPhase(LArPandoraHelper::IsDualPhase());

    if (nWirePlanes != planeSet.size())
      throw cet::exception("LArPandora")
        << " LArPandoraGeometry::LoadGeometry --- geometry description for wire plane(s) missing ";

    if (isDualPhase && (!planeSet.count(geo::kW) || !planeSet.count(geo::kY)))
      throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadGeometry --- dual phase "
                                            "scenario; expect to find w and y views ";

    if (!isDualPhase && (!planeSet.count(geo::kU) || !planeSet.count(geo::kV) ||
                         (planeSet.count(geo::kW) && planeSet.count(geo::kY))))
      throw cet::exception("LArPandora")
        << " LArPandoraGeometry::LoadGeometry --- single phase scenatio; expect to find u and v "
           "views; if there is one further view, it must be w or y ";

    const bool useYPlane((nWirePlanes > 2) && planeSet.count(geo::kY));

    // ATTN: In the dual phase mode, map the wire planes as follows W->U and Y->V.  This mapping was chosen so that the dual phase wire
    // planes, which are inherently induction only, are mapped to induction planes in the single phase geometry.
    const float wirePitchU(theGeometry->WirePitch((isDualPhase ? geo::kW : geo::kU)));
    const float wirePitchV(theGeometry->WirePitch((isDualPhase ? geo::kY : geo::kV)));
    const float wirePitchW((nWirePlanes < 3) ? 0.5f * (wirePitchU + wirePitchV) :
                                               (useYPlane) ? theGeometry->WirePitch(geo::kY) :
                                                             theGeometry->WirePitch(geo::kW));

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

        // ATTN: In dual phase scenario propagate the W->U and Y->V mapping and set wire angle for remaining view to epsilon to
        // avoid identical wire angles clashes (dual phase W and Y wires are horizontal and vertical).  Inside LArSoft the
        // WireAngleToVertical function returns the wire angle to the positive Z axis, but Pandora expects to receive the wire
        // angle to the vertical, hence the conversion.  The fabs() in wireAngleW for the kY case is due to the ICARUS geometry
        // having a wire angle of PI instead of 0.
        const geo::View_t targetViewU(isDualPhase ? geo::kW : geo::kU);
        const geo::View_t targetViewV(isDualPhase ? geo::kY : geo::kV);
        const float wireAngleU(0.5f * M_PI -
                               theGeometry->WireAngleToVertical(targetViewU, itpc1, icstat));
        const float wireAngleV(0.5f * M_PI -
                               theGeometry->WireAngleToVertical(targetViewV, itpc1, icstat));
        const float wireAngleW(
          (nWirePlanes < 3) ?
            std::numeric_limits<float>::epsilon() :
            (useYPlane) ?
            (std::fabs(0.5f * M_PI - theGeometry->WireAngleToVertical(geo::kY, itpc1, icstat))) :
            (0.5f * M_PI - theGeometry->WireAngleToVertical(geo::kW, itpc1, icstat)));

        double localCoord1[3] = {0., 0., 0.};
        double worldCoord1[3] = {0., 0., 0.};
        theTpc1.LocalToWorld(localCoord1, worldCoord1);

        const double min1(worldCoord1[0] - 0.5 * theTpc1.ActiveHalfWidth());
        const double max1(worldCoord1[0] + 0.5 * theTpc1.ActiveHalfWidth());

        float driftMinX(worldCoord1[0] - theTpc1.ActiveHalfWidth());
        float driftMaxX(worldCoord1[0] + theTpc1.ActiveHalfWidth());
        float driftMinY(worldCoord1[1] - theTpc1.ActiveHalfHeight());
        float driftMaxY(worldCoord1[1] + theTpc1.ActiveHalfHeight());
        float driftMinZ(worldCoord1[2] - 0.5f * theTpc1.ActiveLength());
        float driftMaxZ(worldCoord1[2] + 0.5f * theTpc1.ActiveLength());

        const bool isPositiveDrift(theTpc1.DriftDirection() == geo::kPosX);

        UIntSet tpcList;
        tpcList.insert(itpc1);

        // Now identify the other TPCs associated with this drift volume
        for (unsigned int itpc2 = itpc1 + 1; itpc2 < theGeometry->NTPC(icstat); ++itpc2) {
          if (cstatList.end() != cstatList.find(itpc2)) continue;

          const geo::TPCGeo& theTpc2(theGeometry->TPC(itpc2, icstat));

          if (theTpc1.DriftDirection() != theTpc2.DriftDirection()) continue;

          // ATTN: In dual phase scenario propagate the W->U and Y->V mapping as described above.
          const geo::View_t pandoraUView(isDualPhase ? geo::kW : geo::kU);
          const geo::View_t pandoraVView(isDualPhase ? geo::kY : geo::kV);
          const float dThetaU(theGeometry->WireAngleToVertical(pandoraUView, itpc1, icstat) -
                              theGeometry->WireAngleToVertical(pandoraUView, itpc2, icstat));
          const float dThetaV(theGeometry->WireAngleToVertical(pandoraVView, itpc1, icstat) -
                              theGeometry->WireAngleToVertical(pandoraVView, itpc2, icstat));
          const float dThetaW((nWirePlanes < 3) ?
                                std::numeric_limits<float>::epsilon() :
                                (useYPlane) ?
                                (theGeometry->WireAngleToVertical(geo::kY, itpc1, icstat) -
                                 theGeometry->WireAngleToVertical(geo::kY, itpc2, icstat)) :
                                (theGeometry->WireAngleToVertical(geo::kW, itpc1, icstat) -
                                 theGeometry->WireAngleToVertical(geo::kW, itpc2, icstat)));

          if (dThetaU > maxDeltaTheta || dThetaV > maxDeltaTheta || dThetaW > maxDeltaTheta)
            continue;

          double localCoord2[3] = {0., 0., 0.};
          double worldCoord2[3] = {0., 0., 0.};
          theTpc2.LocalToWorld(localCoord2, worldCoord2);

          const double min2(worldCoord2[0] - 0.5 * theTpc2.ActiveHalfWidth());
          const double max2(worldCoord2[0] + 0.5 * theTpc2.ActiveHalfWidth());

          if ((min2 > max1) || (min1 > max2)) continue;

          cstatList.insert(itpc2);
          tpcList.insert(itpc2);

          const float driftMinX2(worldCoord2[0] - theTpc2.ActiveHalfWidth());
          const float driftMaxX2(worldCoord2[0] + theTpc2.ActiveHalfWidth());
          const float driftMinY2(worldCoord2[1] - theTpc2.ActiveHalfHeight());
          const float driftMaxY2(worldCoord2[1] + theTpc2.ActiveHalfHeight());
          const float driftMinZ2(worldCoord2[2] - 0.5f * theTpc2.ActiveLength());
          const float driftMaxZ2(worldCoord2[2] + 0.5f * theTpc2.ActiveLength());

          driftMinX = std::min(driftMinX, driftMinX2);
          driftMaxX = std::max(driftMaxX, driftMaxX2);
          driftMinY = std::min(driftMinY, driftMinY2);
          driftMaxY = std::max(driftMaxY, driftMaxY2);
          driftMinZ = std::min(driftMinZ, driftMinZ2);
          driftMaxZ = std::max(driftMaxZ, driftMaxZ2);
        }

        // Collate the tpc volumes in this drift volume
        LArDaughterDriftVolumeList tpcVolumeList;

        for (const unsigned int itpc : tpcList) {
          tpcVolumeList.emplace_back(icstat, itpc);
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

    // Create daughter drift volumes
    for (const LArDriftVolume& driftVolume : driftVolumeList) {
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

  //------------------------------------------------------------------------------------------------------------------------------------------

  const LArDaughterDriftVolumeList&
  LArDriftVolume::GetTpcVolumeList() const
  {
    return m_tpcVolumeList;
  }

} // namespace lar_pandora
