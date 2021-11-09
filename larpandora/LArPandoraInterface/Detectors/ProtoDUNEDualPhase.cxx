/**
 *  @file   larpandora/LArPandoraInterface/Detectors/ProtoDUNEDualPhase.cxx
 *
 *  @brief  Implementation of the ProtoDUNE dual phase interface
 *
 *  $Log: $
 */

#include "larpandora/LArPandoraInterface/Detectors/ProtoDUNEDualPhase.h"

namespace lar_pandora {

  void
  ProtoDUNEDualPhase::LoadDaughterDetectorGaps(const LArDriftVolume& driftVolume,
                                               const float maxDisplacement,
                                               LArDetectorGapList& listOfGaps) const
  {
    for (LArDaughterDriftVolumeList::const_iterator
           iterDghtr1 = driftVolume.GetTpcVolumeList().begin(),
           iterDghtrEnd1 = driftVolume.GetTpcVolumeList().end();
         iterDghtr1 != iterDghtrEnd1;
         ++iterDghtr1) {
      const LArDaughterDriftVolume& tpcVolume1(*iterDghtr1);

      for (LArDaughterDriftVolumeList::const_iterator
             iterDghtr2 = iterDghtr1,
             iterDghtrEnd2 = driftVolume.GetTpcVolumeList().end();
           iterDghtr2 != iterDghtrEnd2;
           ++iterDghtr2) {
        const LArDaughterDriftVolume& tpcVolume2(*iterDghtr2);

        if (tpcVolume1.GetTpc() == tpcVolume2.GetTpc()) continue;

        const float deltaX(std::fabs(tpcVolume1.GetCenterX() - tpcVolume2.GetCenterX()));
        const float deltaY(std::fabs(tpcVolume1.GetCenterY() - tpcVolume2.GetCenterY()));
        const float deltaZ(std::fabs(tpcVolume1.GetCenterZ() - tpcVolume2.GetCenterZ()));

        const float widthX(0.5f * (tpcVolume1.GetWidthX() + tpcVolume2.GetWidthX()));
        const float widthY(0.5f * (tpcVolume1.GetWidthY() + tpcVolume2.GetWidthY()));
        const float widthZ(0.5f * (tpcVolume1.GetWidthZ() + tpcVolume2.GetWidthZ()));

        const float gapX(deltaX - widthX);
        const float gapY(deltaY - widthY);
        const float gapZ(deltaZ - widthZ);

        const float X1((tpcVolume1.GetCenterX() < tpcVolume2.GetCenterX()) ?
                         (tpcVolume1.GetCenterX() + 0.5f * tpcVolume1.GetWidthX()) :
                         (tpcVolume2.GetCenterX() + 0.5f * tpcVolume2.GetWidthX()));
        const float X2((tpcVolume1.GetCenterX() > tpcVolume2.GetCenterX()) ?
                         (tpcVolume1.GetCenterX() - 0.5f * tpcVolume1.GetWidthX()) :
                         (tpcVolume2.GetCenterX() - 0.5f * tpcVolume2.GetWidthX()));
        const float Y1(std::min((tpcVolume1.GetCenterY() - 0.5f * tpcVolume1.GetWidthY()),
                                (tpcVolume2.GetCenterY() - 0.5f * tpcVolume2.GetWidthY())));
        const float Y2(std::max((tpcVolume1.GetCenterY() + 0.5f * tpcVolume1.GetWidthY()),
                                (tpcVolume2.GetCenterY() + 0.5f * tpcVolume2.GetWidthY())));
        const float Z1(std::min((tpcVolume1.GetCenterZ() - 0.5f * tpcVolume1.GetWidthZ()),
                                (tpcVolume2.GetCenterZ() - 0.5f * tpcVolume2.GetWidthZ())));
        const float Z2(std::max((tpcVolume1.GetCenterZ() + 0.5f * tpcVolume1.GetWidthZ()),
                                (tpcVolume2.GetCenterZ() + 0.5f * tpcVolume2.GetWidthZ())));

        geo::Vector_t gaps(gapX, gapY, gapZ), deltas(deltaX, deltaY, deltaZ);
        if (this->CheckDetectorGapSize(gaps, deltas, maxDisplacement)) {
          geo::Point_t point1(X1, Y1, Z1), point2(X2, Y2, Z2);
          geo::Vector_t widths(widthX, widthY, widthZ);
          listOfGaps.emplace_back(this->CreateDetectorGap(point1, point2, widths));
        }
      }
    }

    return;
  }

} // namespace lar_pandora
