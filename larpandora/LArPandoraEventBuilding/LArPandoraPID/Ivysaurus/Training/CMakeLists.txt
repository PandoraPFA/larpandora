cet_build_plugin(IvysaurusTrainingFiles art::EDAnalyzer
  LIBRARIES
    PRIVATE
      larpandora::LArPandoraUtils
      larpandora::LArPandoraEventBuilding_LArPandoraPID_Ivysaurus_Utils
      larpandora::LArPandoraEventBuilding_LArPandoraPID_Ivysaurus_Managers    
      larpandora::LArPandoraInterface
      larreco::Calorimetry
      larreco::RecoAlg
      larsim::MCCheater_BackTrackerService_service           
      larsim::MCCheater_ParticleInventoryService_service
      larsim::Utils      
)

install_headers()
install_fhicl()
install_source()
