cet_make_library(
  SOURCE 
    IvysaurusUtils.cc 
  LIBRARIES
    PRIVATE
      larcorealg::Geometry
      lardataobj::RecoBase
      larpandora::GetDetectorType
      larpandora::LArPandoraInterface
      larpandora::LArPandoraUtils
      larsim::Simulation 
      larsim::Utils      
)

install_headers()
install_fhicl()
install_source()

