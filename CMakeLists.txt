# ======================================================================
#  larpandora main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larpandora/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================


# use cmake 3.14 or later
cmake_minimum_required (VERSION 3.12...3.14 FATAL_ERROR)

project(larpandora LANGUAGES CXX)

# cetbuildtools contains our cmake modules
SET ( CETBUILDTOOLS_VERSION $ENV{CETBUILDTOOLS_VERSION} )
IF (NOT CETBUILDTOOLS_VERSION)
    MESSAGE (FATAL_ERROR "ERROR: setup cetbuildtools to get the cmake modules")
ENDIF()

set(CMAKE_MODULE_PATH $ENV{CANVAS_ROOT_IO_DIR}/Modules
          $ENV{CETBUILDTOOLS_DIR}/Modules
          ${CMAKE_MODULE_PATH})

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
  )

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_ups_product( larcoreobj )
find_ups_product( larcorealg )
find_ups_product( larcore )
find_ups_product( lardataobj )
find_ups_product( lardataalg )
find_ups_product( lardata )
find_ups_product( larevt )
find_ups_product( larsim )
find_ups_product( larreco )
find_ups_product( nug4)
find_ups_product( art )
find_ups_product(art_root_io)
find_ups_product( pandora )
find_ups_product( larpandoracontent )
find_ups_product( cetbuildtools )
find_ups_product( postgresql )
find_ups_product( eigen )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

cet_find_library( PANDORASDK NAMES PandoraSDK PATHS ENV PANDORA_LIB )
cet_find_library( PANDORAMONITORING NAMES PandoraMonitoring PATHS ENV PANDORA_LIB )

# find larpandoracontent headers if building at the same time
#message(STATUS "larpandora: checking for MRB_SOURCE")
set( mrb_source $ENV{MRB_SOURCE} )
if( LARPANDORACONTENT_INC )
  #message(STATUS "larpandora: found LARPANDORACONTENT_INC: $ENV{LARPANDORACONTENT_INC}")
  # note that mrb v1_07_00 or later will define LARPANDORACONTENT_INC
  if( mrb_source )
    include_directories( ${LARPANDORACONTENT_INC} )
  endif()
else()
  if( mrb_source )
    include_directories( ${mrb_source}/larpandoracontent )
  else()
    message(FATAL_ERROR "LARPANDORACONTENT_INC is undefined, please setup larpandoracontent")
  endif()
endif()

# source
add_subdirectory(larpandora)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)

