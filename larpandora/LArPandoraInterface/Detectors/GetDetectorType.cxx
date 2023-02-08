/**
 *  @file   larpandora/LArPandoraInterface/Detectors/LArPandoraDetectorType.cxx
 *
 *  @brief  Implementation of the interface for handling detector-specific details, as well as some helper functions
 *
 *  $Log: $
 */

#include "larpandora/LArPandoraInterface/Detectors/GetDetectorType.h"
#include "larpandora/LArPandoraInterface/Detectors/DUNEFarDetVDThreeView.h"
#include "larpandora/LArPandoraInterface/Detectors/DUNEFarDetVDThreeView30DegDriftY.h"
#include "larpandora/LArPandoraInterface/Detectors/ICARUS.h"
#include "larpandora/LArPandoraInterface/Detectors/LArPandoraDetectorType.h"
#include "larpandora/LArPandoraInterface/Detectors/ProtoDUNEDualPhase.h"
#include "larpandora/LArPandoraInterface/Detectors/VintageLArTPCThreeView.h"

#include "cetlib_except/exception.h"

#include <limits>
#include <set>

namespace lar_pandora {

  LArPandoraDetectorType* detector_functions::GetDetectorType()
  {
    art::ServiceHandle<geo::Geometry const> geo;

    std::set<short int> driftDirectionSet;
    for (geo::CryostatGeo const& cryo: geo->Iterate<geo::CryostatGeo>()){
        for (geo::TPCGeo const& TPC: cryo.IterateTPCs()) {
            (void)driftDirectionSet.insert(TPC.DetectDriftDirection());
        }
    }
    const int driftDirectionCount( (driftDirectionSet.count(1) || driftDirectionSet.count(-1)) +
                                   (driftDirectionSet.count(2) || driftDirectionSet.count(-2)) +
                                   (driftDirectionSet.count(3) || driftDirectionSet.count(-3)) +
                                   driftDirectionSet.count(0) );

    if (driftDirectionCount > 1)
        throw cet::exception("LArPandora") << "LArPandoraDetectorType::GetDetectorType -- unable to "
                                              "determine drift direction due to erroneous count.  "
                                              "Present drift directions: \n"
                                              "-- +x: " <<  driftDirectionSet.count(1)  << "\n"
                                              "-- -x: " <<  driftDirectionSet.count(-1) << "\n"
                                              "-- +y: " <<  driftDirectionSet.count(2)  << "\n"
                                              "-- -y: " <<  driftDirectionSet.count(-2) << "\n"
                                              "-- +z: " <<  driftDirectionSet.count(3)  << "\n"
                                              "-- -z: " <<  driftDirectionSet.count(-3) << "\n"
                                              "-- ??: " <<  driftDirectionSet.count(0)  << "\n";

    const unsigned int nPlanes(geo->MaxPlanes());
    std::set<geo::_plane_proj> planeSet;
    for (unsigned int iPlane = 0; iPlane < nPlanes; ++iPlane)
      (void)planeSet.insert(geo->TPC().Plane(iPlane).View());

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

} // namespace lar_pandora
