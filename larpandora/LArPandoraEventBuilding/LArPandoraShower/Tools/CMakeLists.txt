cet_make_library(LIBRARY_NAME ShowerTool INTERFACE
  SOURCE IShowerTool.h
  LIBRARIES INTERFACE
  larpandora::LArPandoraEventBuilding_LArPandoraShower_Algs
  larpandora::ShowerProducedPtrs
  larpandora::ShowerElements
)

cet_write_plugin_builder(larpandora::ShowerTool art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  larpandora::ShowerTool
)

include(larpandora::ShowerTool)

cet_build_plugin(Shower2DLinearRegressionTrackHitFinder larpandora::ShowerTool
  LIBRARIES PRIVATE
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardataobj::RecoBase
)

cet_build_plugin(Shower3DCylinderTrackHitFinder larpandora::ShowerTool
  LIBRARIES PRIVATE
)

cet_build_plugin(ShowerBayesianTrucatingdEdx larpandora::ShowerTool
  LIBRARIES PRIVATE
  ROOT::Hist
  ROOT::RIO
)

cet_build_plugin(ShowerDirectionTopologyDecision larpandora::ShowerTool
  LIBRARIES PRIVATE
)

cet_build_plugin(ShowerExampleTool larpandora::ShowerTool
  LIBRARIES PRIVATE
  lardataobj::RecoBase
)

cet_build_plugin(ShowerIncrementalTrackHitFinder larpandora::ShowerTool
  LIBRARIES PRIVATE
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardataalg::DetectorInfo
  lardataobj::RecoBase
  ROOT::Gpad
  ROOT::Hist
  ROOT::RIO
)

cet_build_plugin(ShowerLengthPercentile larpandora::ShowerTool
  LIBRARIES PRIVATE
)

cet_build_plugin(ShowerLinearEnergy larpandora::ShowerTool
  LIBRARIES PRIVATE
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardataobj::RecoBase
)

cet_build_plugin(ShowerNumElectronsEnergy larpandora::ShowerTool
  LIBRARIES PRIVATE
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larreco::Calorimetry
  lardataobj::RecoBase
  larcoreobj::headers
)

cet_build_plugin(ShowerPCADirection larpandora::ShowerTool
  LIBRARIES PRIVATE
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardataalg::DetectorInfo
  lardataobj::RecoBase
  ROOT::Hist
  Eigen3::Eigen
)

cet_build_plugin(ShowerPCAEigenvalueLength larpandora::ShowerTool
  LIBRARIES PRIVATE
)

cet_build_plugin(ShowerPCAPropergationStartPosition larpandora::ShowerTool
  LIBRARIES PRIVATE
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
)

cet_build_plugin(ShowerPFPVertexStartPosition larpandora::ShowerTool
  LIBRARIES PRIVATE
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardataobj::RecoBase
)

cet_build_plugin(ShowerPandoraSlidingFitTrackFinder larpandora::ShowerTool
  LIBRARIES PRIVATE
  larpandora::GetDetectorType
  larpandora::LArPandoraInterface_Detectors
  larpandoracontent::LArPandoraContent
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardataobj::RecoBase
)

cet_build_plugin(ShowerSkeletonTool larpandora::ShowerTool
  LIBRARIES PRIVATE
)

cet_build_plugin(ShowerTrackColinearTrajPointDirection larpandora::ShowerTool
  LIBRARIES PRIVATE
  lardataobj::RecoBase
)

cet_build_plugin(ShowerTrackDirection larpandora::ShowerTool
  LIBRARIES PRIVATE
  lardataobj::RecoBase
)

cet_build_plugin(ShowerTrackHitDirection larpandora::ShowerTool
  LIBRARIES PRIVATE
  lardataobj::RecoBase
)

cet_build_plugin(ShowerTrackPCADirection larpandora::ShowerTool
  LIBRARIES PRIVATE
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardataalg::DetectorInfo
)

cet_build_plugin(ShowerTrackSpacePointDirection larpandora::ShowerTool
  LIBRARIES PRIVATE
  lardataobj::RecoBase
)

cet_build_plugin(ShowerTrackStartPosition larpandora::ShowerTool
  LIBRARIES PRIVATE
  lardataobj::RecoBase
)

cet_build_plugin(ShowerTrackTrajPointDirection larpandora::ShowerTool
  LIBRARIES PRIVATE
  lardataobj::RecoBase
)

cet_build_plugin(ShowerTrackTrajToSpacePoint larpandora::ShowerTool
  LIBRARIES PRIVATE
  lardataobj::RecoBase
)

cet_build_plugin(ShowerTrajPointdEdx larpandora::ShowerTool
  LIBRARIES PRIVATE
  larreco::Calorimetry
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardataobj::RecoBase
)

cet_build_plugin(ShowerUnidirectiondEdx larpandora::ShowerTool
  LIBRARIES PRIVATE
  larreco::Calorimetry
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larpandora::LArPandoraInterface

)

add_subdirectory(Cheating)

install_headers()
install_source()
install_fhicl()
