include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

add_subdirectory(Evaluator)  
add_subdirectory(Managers)
add_subdirectory(Training)
add_subdirectory(Utils)

install_headers()
install_fhicl()
install_source()

