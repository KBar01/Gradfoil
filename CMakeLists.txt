cmake_minimum_required(VERSION "3.5.0")
project("gradFoil")

# Include my projects header files
include_directories(/src/include)

#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

# Find Eigen (installed in system)
find_package(Eigen3 REQUIRED)
include_directories(/usr/include/eigen3)

# Find CoDiPack
find_package(CoDiPack REQUIRED)
include_directories(/usr/include/codi)

# Find JSON stuff
include_directories(/usr/local/include)

# Gather all source files
file(GLOB_RECURSE SOURCES src/*.cpp)


########## Defining two executables ##################
# --- Standard executable ---
add_executable(CFoil_fwd main.cpp ${SOURCES})
target_link_libraries(CFoil_fwd Eigen3::Eigen)
target_compile_definitions(CFoil_fwd PRIVATE -DCODI_EnableEigen)

# --- Codipack-enabled executable ---
add_executable(CFoil_AD main.cpp ${SOURCES})
target_link_libraries(CFoil_AD Eigen3::Eigen CoDiPack)
target_compile_definitions(CFoil_AD PRIVATE -DUSE_CODIPACK -DCODI_EnableEigen)
