

cmake_minimum_required(VERSION 3.14)
project(gradFoil)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(CMAKE_CXX_FLAGS_RELEASE "-O2 -fno-fast-math")
  #-fno-strict-aliasing
endif()

option(DO_BL_GRADIENT "Enable boundary layer gradient calculation" ON)
option(DO_SOUND "Enable calculation of TE sound" OFF)

if(DO_BL_GRADIENT)
  set(DO_BL_GRADIENT_DEFINE 1)
else()
  set(DO_BL_GRADIENT_DEFINE 0)
endif()

if(DO_SOUND)
  set(DO_SOUND_DEFINE 1)
else()
  set(DO_SOUND_DEFINE 0)
endif()


include(FetchContent)

#############################################
# Fetch Eigen (header-only)
#############################################
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0
)
FetchContent_MakeAvailable(eigen)

#############################################
# Fetch CoDiPack (header-only)
#############################################

FetchContent_Declare(
  codipack
  GIT_REPOSITORY https://github.com/SciCompKL/CoDiPack.git
  GIT_TAG master  # use a specific commit/tag if needed
)
FetchContent_MakeAvailable(codipack)

#############################################
# Fetch json (header-only)
#############################################

FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz)
FetchContent_MakeAvailable(json)

#############################################
# Include your own headers
#############################################
include_directories(${PROJECT_SOURCE_DIR}/src/include)

# Include CoDiPack manually (header-only)
#include_directories(${codipack_SOURCE_DIR}/include)

#############################################
# Gather source files
#############################################
file(GLOB_RECURSE SOURCES src/*.cpp)
list(REMOVE_ITEM SOURCES "${CMAKE_SOURCE_DIR}/src/main.cpp")
#############################################
# Build Executables
#############################################
add_executable(CFoil_fwd src/main.cpp ${SOURCES})
target_link_libraries(CFoil_fwd Eigen3::Eigen CoDiPack nlohmann_json::nlohmann_json)
target_compile_definitions(CFoil_fwd PRIVATE
  -DCODI_EnableEigen
  -DDO_BL_GRADIENT=${DO_BL_GRADIENT_DEFINE}
  -DDO_SOUND=${DO_SOUND_DEFINE}
)

add_executable(CFoil_AD src/main.cpp ${SOURCES})
target_link_libraries(CFoil_AD PRIVATE Eigen3::Eigen CoDiPack nlohmann_json::nlohmann_json)

target_include_directories(CFoil_AD PRIVATE
    ${codipack_SOURCE_DIR}/include
)

target_compile_definitions(CFoil_AD PRIVATE
  -DUSE_CODIPACK
  -DCODI_EnableEigen 
  -DCODI_ForcedInlines
  -DDO_BL_GRADIENT=${DO_BL_GRADIENT_DEFINE}
  -DDO_SOUND=${DO_SOUND_DEFINE}
)