find_ups_product(tensorflow)
find_ups_product(protobuf)

cet_add_compiler_flags(CXX -Wno-pedantic)

cet_make_library(
  SOURCE 
    IvysaurusEvaluator.cxx
  LIBRARIES
    PRIVATE
      $ENV{TENSORFLOW_LIB}/libtensorflow_cc.so
      $ENV{TENSORFLOW_LIB}/libtensorflow_framework.so
      $ENV{PROTOBUF_LIB}/libprotobuf.so
      pthread
      Eigen3::Eigen
      larcorealg::Geometry
      larcore::Geometry_Geometry_service
      larsim::Simulation 
      lardataobj::Simulation
      lardataobj::RawData
      lardataobj::RecoBase
      lardataobj::AnalysisBase
      larreco::Calorimetry
      lardata::Utilities
      lardataalg::DetectorInfo
      nusimdata::SimulationBase
      larsim::MCCheater_BackTrackerService_service           
      larsim::MCCheater_ParticleInventoryService_service
      ROOT::Core
      art::Persistency_Common
      art::Persistency_Provenance
      art::Utilities
      canvas::canvas
      messagefacility::MF_MessageLogger
      fhiclcpp::fhiclcpp
      cetlib::cetlib
      larpandora::LArPandoraInterface
      larreco::RecoAlg
      larpandora::LArPandoraUtils
      larpandora::LArPandoraEventBuilding_LArPandoraPID_Ivysaurus_Managers
)

install_headers()
install_fhicl()
install_source()
