include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

if(( DEFINED ENV{TENSORFLOW_DIR} ) AND ( DEFINED ENV{TRITON_DIR} ))
  add_subdirectory(Managers)
  add_subdirectory(Modules)
  add_subdirectory(TensorFlow)  
  add_subdirectory(Utils)
endif ()

install_headers()
install_fhicl()
install_source()

