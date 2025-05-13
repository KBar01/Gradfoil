

cmake_minimum_required(VERSION 3.14)
project(gradFoil)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)

# Output binaries to source directory (for setup.py packaging)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

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
target_link_libraries(CFoil_fwd Eigen3::Eigen)
#target_compile_definitions(CFoil_fwd PRIVATE -DCODI_EnableEigen)

add_executable(CFoil_AD src/main.cpp ${SOURCES})
target_link_libraries(CFoil_AD PRIVATE Eigen3::Eigen)

target_include_directories(CFoil_AD PRIVATE
    ${codipack_SOURCE_DIR}/include
)

target_compile_definitions(CFoil_AD PRIVATE -DUSE_CODIPACK -DCODI_EnableEigen -DCODI_ForcedInlines)