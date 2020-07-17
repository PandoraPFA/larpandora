include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{LARPANDORACONTENT_INC} )

add_subdirectory(Tools)
add_subdirectory(Algs)

art_make(
  MODULE_LIBRARIES
  larpandora_LArPandoraEventBuilding
  larpandora_LArPandoraEventBuilding_LArPandoraShower_Algs
  )


install_fhicl()
install_source()
