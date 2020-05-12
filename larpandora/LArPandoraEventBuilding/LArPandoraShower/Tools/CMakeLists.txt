include_directories( $ENV{PANDORA_INC} )

add_subdirectory(Cheating)

art_make(TOOL_LIBRARIES
  larpandora_LArPandoraEventBuilding_LArPandoraShower_Algs
  larreco_Calorimetry
  LArPandoraContent
  )

install_headers()
install_source()
install_fhicl()
