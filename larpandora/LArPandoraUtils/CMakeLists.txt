cet_make_library(
  LIBRARY_NAME PandoraUtilsBase INTERFACE
  SOURCE PandoraUtilsBase.h
  LIBRARIES INTERFACE
    canvas::canvas  
)

cet_make_library(
  SOURCE 
    PandoraClusterUtils.cxx
    PandoraEventUtils.cxx
    PandoraHitUtils.cxx        
    PandoraPFParticleUtils.cxx
    PandoraShowerUtils.cxx
    PandoraSliceUtils.cxx
    PandoraSpacePointUtils.cxx    
    PandoraTrackUtils.cxx    
  LIBRARIES
    PRIVATE
      art::Persistency_Common
      art::Persistency_Provenance    
      lardata::Utilities
      lardataobj::Simulation  
      lardataobj::RecoBase
      larpandora::LArPandoraInterface
      larpandora::PandoraUtilsBase
)

install_headers()
install_fhicl()
install_source()

