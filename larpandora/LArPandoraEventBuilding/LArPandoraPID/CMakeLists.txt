add_subdirectory(LArPandoraPIDData)
add_subdirectory(Ivysaurus)

find_ups_product(protobuf)

cet_add_compiler_flags(CXX -Wno-pedantic)

cet_build_plugin(PandoraPID art::EDProducer
  LIBRARIES
    PRIVATE
      art::Utilities
      Eigen3::Eigen    
      lardataobj::RecoBase
      larpandora::LArPandoraUtils
      larpandora::LArPandoraEventBuilding_LArPandoraPID_Ivysaurus_Managers    
      larpandora::LArPandoraEventBuilding_LArPandoraPID_Ivysaurus_Evaluator
      larpandora::LArPandoraEventBuilding_LArPandoraPID_LArPandoraPIDData      
      larpandora::LArPandoraInterface
      larreco::Calorimetry
      torch
      torch_cpu   
)


install_fhicl()
install_source()

