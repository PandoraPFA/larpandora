include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{LARPANDORACONTENT_INC} )

add_subdirectory(Tools)
add_subdirectory(Algs)

art_make(
  LIB_LIBRARIES larcorealg_Geometry
  larcore_Geometry_Geometry_service
  larsim_Simulation lardataobj_Simulation
  lardataobj_RawData
  lardataobj_RecoBase
  lardataobj_AnalysisBase
  lardata_ArtDataHelper
  lardata_Utilities
  larreco_RecoAlg_ClusterRecoUtil
  larreco_RecoAlg
  larreco_ClusterFinder
  larreco_Calorimetry
  ${PANDORASDK}
  ${PANDORAMONITORING}
  LArPandoraContent
  larpandora_LArPandoraInterface
  nusimdata_SimulationBase
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_BASIC}
  ${ART_ROOT_IO_TFILE_SUPPORT} ROOT::Core
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  art_Persistency_Common
  art_Persistency_Provenance
  art_Utilities
  canvas
  ${MF_MESSAGELOGGER}
  ${FHICLCPP}
  cetlib cetlib_except
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ROOT::Geom
  ${ROOT_BASIC_LIB_LIST}
  ROOT::GenVector
  MODULE_LIBRARIES larpandora_LArPandoraEventBuilding
  larpandora_LArPandoraEventBuilding_LArPandoraShower_Algs
  )


install_headers()
install_fhicl()
install_source()
