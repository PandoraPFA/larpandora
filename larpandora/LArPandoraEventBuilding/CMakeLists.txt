cet_make_library(LIBRARY_NAME Slice INTERFACE
  SOURCE Slice.h
  LIBRARIES INTERFACE
  larpandora::LArPandoraInterface 
)

cet_make_library(LIBRARY_NAME SliceIdTool INTERFACE
  SOURCE SliceIdBaseTool.h
  LIBRARIES INTERFACE
  larpandora::Slice
)

cet_write_plugin_builder(larpandora::SliceIdTool art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  larpandora::SliceIdTool
)

cet_make_library(SOURCE
  LArPandoraEvent.cxx
  LArPandoraSliceIdHelper.cxx
  LIBRARIES
  PUBLIC
  larpandora::LArPandoraInterface
  lardataobj::AnalysisBase
  lardataobj::RecoBase
  art::Framework_Principal
  art::Persistency_Common
  canvas::canvas
  PRIVATE
  nusimdata::SimulationBase
)

add_subdirectory(LArPandoraPID)
add_subdirectory(LArPandoraShower)
add_subdirectory(SliceIdTools)

cet_build_plugin(CollectionSplitting art::EDProducer
  LIBRARIES PRIVATE
  larpandora::LArPandoraEventBuilding
  lardataobj::RecoBase
  art::Framework_Principal
  fhiclcpp::fhiclcpp
)

cet_build_plugin(LArPandoraEventDump art::EDAnalyzer
  LIBRARIES PRIVATE
  larpandora::LArPandoraInterface
  lardata::AssociationUtil
  lardataobj::RecoBase
  art::Framework_Principal
  fhiclcpp::fhiclcpp
)

cet_build_plugin(LArPandoraExternalEventBuilding art::EDProducer
  LIBRARIES PRIVATE
  larpandora::LArPandoraEventBuilding
  larpandora::Slice
  larpandora::SliceIdTool
  larpandora::LArPandoraInterface
  lardata::AssociationUtil
  lardataobj::RecoBase
  art::Framework_Principal
  art_plugin_support::toolMaker
  canvas::canvas
  fhiclcpp::fhiclcpp
  PandoraPFA::PandoraSDK
)

cet_build_plugin(LArPandoraTrackCreation art::EDProducer
  LIBRARIES PRIVATE
  larpandora::GetDetectorType
  larpandora::LArPandoraInterface_Detectors
  larpandora::LArPandoraInterface
  larpandoracontent::LArPandoraContent
  lardata::AssociationUtil
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art::Framework_Principal
  art::Persistency_Common
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
)

add_subdirectory(fcl)

install_headers()
install_fhicl()
install_source()
