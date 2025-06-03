

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
# Common include directory for all targets (if needed)
set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/src/include)
set(NOISE_INCLUDES ${PROJECT_SOURCE_DIR}/src/noise_includes)
#include_directories(${PROJECT_SOURCE_DIR}/src/include)

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


# Issue with custom derivative rules and codipack vector mode means
# generating 2 seperate AD scripts when doing OASPL gradient
#
# When DO_SOUND= true : 
# Cfoil_Noise returns just gradient of OASPL and Cfoil_AD retuns
# gradient of CL and CD



################## FWD MODE: includes everything ###############################################
add_executable(CFoil_fwd src/main.cpp ${SOURCES})

target_link_libraries(CFoil_fwd Eigen3::Eigen CoDiPack nlohmann_json::nlohmann_json)

target_include_directories(CFoil_fwd PRIVATE
 ${COMMON_INCLUDES}
 ${NOISE_INCLUDES}
)

target_compile_definitions(CFoil_fwd PRIVATE
  -DCODI_EnableEigen
  -DDO_BL_GRADIENT=${DO_BL_GRADIENT_DEFINE}
  -DDO_SOUND=${DO_SOUND_DEFINE}
)


################### NOISE AD CODE (only if DO_SOUND=ON) ##########################################
if(DO_SOUND)
  
  add_executable(CFoil_Noise src/main.cpp ${SOURCES})
  
  target_link_libraries(CFoil_Noise PRIVATE Eigen3::Eigen CoDiPack nlohmann_json::nlohmann_json)
  
  target_include_directories(CFoil_Noise PRIVATE
  
    ${COMMON_INCLUDES}
    ${NOISE_INCLUDES} 
    ${codipack_SOURCE_DIR}/include
  )
  
  target_compile_definitions(CFoil_Noise PRIVATE
    -DUSE_CODIPACK
    -DCODI_EnableEigen
    -DCODI_ForcedInlines
    -DDO_BL_GRADIENT=0
    -DDO_SOUND=1
  )
endif()


############## Main AD CODE, retuns CL,CD unless BL_GRADIENT=ON, return BL states grads ##########
list(REMOVE_ITEM SOURCES
 "${CMAKE_SOURCE_DIR}/src/amiet.cpp"
 "${CMAKE_SOURCE_DIR}/src/error_package.cpp"
 "${CMAKE_SOURCE_DIR}/src/fresenel.cpp"
 "${CMAKE_SOURCE_DIR}/src/sound.cpp"
 "${CMAKE_SOURCE_DIR}/src/WPSmodels.cpp")

if (DO_SOUND)
  
  add_executable(CFoil_AD src/main.cpp ${SOURCES})
  
  target_link_libraries(CFoil_AD PRIVATE Eigen3::Eigen CoDiPack nlohmann_json::nlohmann_json)

  target_include_directories(CFoil_AD PRIVATE
    ${COMMON_INCLUDES}
    ${codipack_SOURCE_DIR}/include
  )

  target_compile_definitions(CFoil_AD PRIVATE
    -DUSE_CODIPACK
    -DCODI_EnableEigen 
    -DCODI_ForcedInlines
    -DDO_BL_GRADIENT=0
    -DDO_SOUND=0
  )
else()
  
  add_executable(CFoil_AD src/main.cpp ${SOURCES})
  
  target_link_libraries(CFoil_AD PRIVATE Eigen3::Eigen CoDiPack nlohmann_json::nlohmann_json)

  target_include_directories(CFoil_AD PRIVATE
    ${COMMON_INCLUDES}
    ${codipack_SOURCE_DIR}/include
  )

  target_compile_definitions(CFoil_AD PRIVATE
    -DUSE_CODIPACK
    -DCODI_EnableEigen 
    -DCODI_ForcedInlines
    -DDO_BL_GRADIENT=${DO_BL_GRADIENT_DEFINE}
    -DDO_SOUND=0
  )
endif()