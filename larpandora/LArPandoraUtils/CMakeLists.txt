cet_make_library(
  LIBRARY_NAME PandoraUtilsBase INTERFACE
  SOURCE PandoraUtilsBase.h
  LIBRARIES INTERFACE
    canvas::canvas  
)

cet_make_library(
  SOURCE 
    PandoraEventUtils.cxx
    PandoraClusterUtils.cxx
    PandoraHitUtils.cxx        
    PandoraPFParticleUtils.cxx
    PandoraShowerUtils.cxx
    PandoraSliceUtils.cxx
    PandoraSpacePointUtils.cxx    
    PandoraTrackUtils.cxx    
  LIBRARIES
    PRIVATE
      larpandora::PandoraUtilsBase
      larcorealg::Geometry
      larcore::Geometry_Geometry_service
      larsim::Simulation lardataobj::Simulation
      lardataalg::DetectorInfo
      lardataobj::RawData
      lardataobj::RecoBase
      lardataobj::AnalysisBase
      lardata::Utilities
      larpandora::LArPandoraInterface
      nusimdata::SimulationBase
      ROOT::Core
      art::Persistency_Common
      art::Persistency_Provenance
      art::Utilities
      canvas::canvas
      messagefacility::MF_MessageLogger
      fhiclcpp::fhiclcpp
      cetlib::cetlib          
)



install_headers()
install_fhicl()
install_source()

