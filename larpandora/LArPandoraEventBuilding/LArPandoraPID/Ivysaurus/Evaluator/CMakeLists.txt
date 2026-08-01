cet_make_library(
  SOURCE 
    IvysaurusEvaluator.cxx
  LIBRARIES
    PRIVATE
      torch
      torch_cpu
      Eigen3::Eigen
      lardataobj::RecoBase
      larreco::RecoAlg
      larpandora::LArPandoraUtils
      larpandora::LArPandoraEventBuilding_LArPandoraPID_Ivysaurus_Managers
)

install_headers()
install_fhicl()
install_source()
