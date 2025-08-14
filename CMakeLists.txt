# ======================================================================
#  larpandora main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larpandora/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

find_package(cetmodules REQUIRED)
project(larpandora VERSION 10.00.25 LANGUAGES CXX)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
  )

cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

find_package(art REQUIRED EXPORT)
find_package(art_root_io REQUIRED EXPORT)
find_package(canvas REQUIRED EXPORT)
find_package(messagefacility REQUIRED EXPORT)
find_package(fhiclcpp REQUIRED EXPORT)
find_package(cetlib REQUIRED EXPORT)
find_package(cetlib_except REQUIRED EXPORT)

find_package(nusimdata REQUIRED EXPORT)

find_package(larcore REQUIRED EXPORT)
find_package(larcorealg REQUIRED)
find_package(larcoreobj REQUIRED EXPORT)
find_package(lardata REQUIRED)
find_package(lardataalg REQUIRED)
find_package(lardataobj REQUIRED EXPORT)
find_package(larevt REQUIRED)
find_package(larpandoracontent REQUIRED EXPORT) # Must precede Pandora*
find_package(larreco REQUIRED)
find_package(larsim REQUIRED EXPORT)

find_package(Eigen3 3.3 REQUIRED)
find_package(PandoraSDK REQUIRED EXPORT)
find_package(ROOT COMPONENTS Core Gpad Graf3d Hist Physics RIO Tree REQUIRED)

find_package(Torch QUIET EXPORT)
if (Torch_FOUND)
  option(LARPANDORA_LIBTORCH "Flag for building with Pandora's LibTorch-aware algorithms" ON)
else()
  option(LARPANDORA_LIBTORCH "Flag for building with Pandora's LibTorch-aware algorithms" OFF)
endif()
if (LARPANDORA_LIBTORCH)
  message(STATUS "Building WITH deep learning support")
else()
  message(STATUS "Building WITHOUT deep learning support")
endif()

# macros for dictionary and simple_plugin
include(CetMake)
include(BuildPlugins)

cet_cmake_module_directories(Modules BINARY)

# source
add_subdirectory(larpandora)

# packaging utility
cet_cmake_config()
