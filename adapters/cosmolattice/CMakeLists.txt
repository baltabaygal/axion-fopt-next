cmake_minimum_required(VERSION 3.20)
project(axion_fopt_cosmolattice LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
set(COSMOLATTICE_SOURCE "" CACHE PATH "Unmodified CosmoLattice v2.0.0 checkout")
if(NOT EXISTS "${COSMOLATTICE_SOURCE}/include/CosmoInterface/abstractmodel.h")
  message(FATAL_ERROR "Set COSMOLATTICE_SOURCE to the official v2.0.0 checkout")
endif()
find_package(Python3 REQUIRED COMPONENTS Interpreter)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/prepare_overlay.py")
execute_process(COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/prepare_overlay.py"
  "${COSMOLATTICE_SOURCE}" "${CMAKE_CURRENT_BINARY_DIR}/overlay"
  COMMAND_ERROR_IS_FATAL ANY)
# One CPU process only. Keep Serial as the default; OpenMP is explicit and
# validated separately before use in long pocket experiments.
option(AFOPT_KOKKOS_OPENMP "Use Kokkos OpenMP in a single CPU process" OFF)
set(MPI OFF CACHE BOOL "" FORCE)
set(HDF5 OFF CACHE BOOL "" FORCE)
set(CUDA OFF CACHE BOOL "" FORCE)
set(HIP OFF CACHE BOOL "" FORCE)
set(OPENMP ${AFOPT_KOKKOS_OPENMP} CACHE BOOL "" FORCE)
set(PTHREADS OFF CACHE BOOL "" FORCE)
if(AFOPT_KOKKOS_OPENMP)
  set(_afopt_serial OFF)
  set(_afopt_backend "OpenMP")
else()
  set(_afopt_serial ON)
  set(_afopt_backend "Serial")
endif()
set(NOTHREADING ${_afopt_serial} CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_THREADS OFF CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_OPENMP ${AFOPT_KOKKOS_OPENMP} CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_SERIAL ${_afopt_serial} CACHE BOOL "" FORCE)
add_subdirectory("${COSMOLATTICE_SOURCE}" upstream EXCLUDE_FROM_ALL)
# Infinity is part of the arrival-map contract. Upstream Release enables
# finite-only fast math, which invalidates infinity comparisons and validation.
get_target_property(_templ_options TempLat INTERFACE_COMPILE_OPTIONS)
list(REMOVE_ITEM _templ_options -ffast-math -Wno-nan-infinity-disabled)
set_property(TARGET TempLat PROPERTY INTERFACE_COMPILE_OPTIONS "${_templ_options}")
FetchContent_GetProperties(TempLat SOURCE_DIR _templ_source)
execute_process(COMMAND git -C "${_templ_source}" rev-parse HEAD
  OUTPUT_VARIABLE _templ_commit OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
if(NOT _templ_commit STREQUAL "f5239eebdfaa159e45b38b89d7ef8ef6e2222f4b")
  message(FATAL_ERROR "Require official TempLat v1.0.0 commit")
endif()
file(GENERATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/build-info.json" CONTENT
  "{\n  \"cosmolattice_commit\": \"0dad16ec2c894a7c6abbea0e98734126b3470b36\",\n  \"templat_commit\": \"${_templ_commit}\",\n  \"kokkos\": \"5.1.1\",\n  \"compiler\": \"${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}\",\n  \"backend\": \"Kokkos ${_afopt_backend}; FFTW; MPI OFF; strict floating point\"\n}\n")
add_executable(axion_fopt_evolve evolve.cpp)
target_include_directories(axion_fopt_evolve BEFORE PRIVATE
  "${CMAKE_CURRENT_BINARY_DIR}/overlay" "${COSMOLATTICE_SOURCE}/include")
target_link_libraries(axion_fopt_evolve PRIVATE TempLat::TempLat)
