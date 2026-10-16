cmake_minimum_required(VERSION 3.20)
project(pockets LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/cosmolattice" adapter)
add_executable(pocket_evolve evolve.cpp)
target_include_directories(pocket_evolve BEFORE PRIVATE
  "${CMAKE_CURRENT_BINARY_DIR}/adapter/overlay" "${COSMOLATTICE_SOURCE}/include")
target_link_libraries(pocket_evolve PRIVATE TempLat::TempLat)
