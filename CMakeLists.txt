cmake_minimum_required(VERSION 3.20)

project(FTDQMC_Studio LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(FTDQMC_ENABLE_MPI "Enable MPI" ON)
option(FTDQMC_ENABLE_PETSC "Enable PETSc" OFF)
option(FTDQMC_ENABLE_SLEPC "Enable SLEPc" OFF)
option(FTDQMC_ENABLE_HDF5 "Enable HDF5" OFF)

# Backend default: prefer Accelerate on macOS, MKL on Linux, OpenBLAS otherwise.
if(NOT DEFINED FTDQMC_LINALG_BACKEND)
  if(APPLE)
    set(FTDQMC_LINALG_BACKEND "ACCELERATE" CACHE STRING "Linear algebra backend (MKL|ACCELERATE|OPENBLAS|PETSC)")
  elseif(UNIX)
    set(FTDQMC_LINALG_BACKEND "MKL" CACHE STRING "Linear algebra backend (MKL|ACCELERATE|OPENBLAS|PETSC)")
  else()
    set(FTDQMC_LINALG_BACKEND "OPENBLAS" CACHE STRING "Linear algebra backend (MKL|ACCELERATE|OPENBLAS|PETSC)")
  endif()
endif()

string(TOUPPER "${FTDQMC_LINALG_BACKEND}" FTDQMC_LINALG_BACKEND_UPPER)
if(FTDQMC_LINALG_BACKEND_UPPER STREQUAL "PETSC" OR FTDQMC_LINALG_BACKEND_UPPER STREQUAL "SLEPC")
  message(FATAL_ERROR "FTDQMC_LINALG_BACKEND cannot be PETSC/SLEPC. Use MKL/ACCELERATE/OPENBLAS.")
endif()

if(FTDQMC_ENABLE_SLEPC AND NOT FTDQMC_ENABLE_PETSC)
  message(FATAL_ERROR "FTDQMC_ENABLE_SLEPC requires FTDQMC_ENABLE_PETSC=ON")
endif()

include(FetchContent)

FetchContent_Declare(
  tomlplusplus
  GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
  GIT_TAG v3.4.0
)
FetchContent_MakeAvailable(tomlplusplus)

if(FTDQMC_ENABLE_MPI)
  find_package(MPI REQUIRED)
endif()

if(FTDQMC_ENABLE_PETSC)
  # Try PETSc CMake config first. If missing, fall back to pkg-config.
  find_package(PETSc QUIET)
  if(NOT PETSc_FOUND)
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
      if(DEFINED ENV{PETSC_DIR} AND DEFINED ENV{PETSC_ARCH})
        set(ENV{PKG_CONFIG_PATH} "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
      endif()
      pkg_check_modules(PETSC QUIET PETSc petsc)
      if(PETSC_FOUND)
        add_library(PETSc::PETSc INTERFACE IMPORTED)
        target_include_directories(PETSc::PETSc INTERFACE ${PETSC_INCLUDE_DIRS})
        target_link_libraries(PETSc::PETSc INTERFACE ${PETSC_LINK_LIBRARIES})
        set(PETSc_FOUND TRUE)
      endif()
    endif()
  endif()
  if(NOT PETSc_FOUND AND NOT PETSC_FOUND)
    message(FATAL_ERROR "PETSc not found. Set PETSc_DIR or PETSC_DIR/PETSC_ARCH for pkg-config fallback.")
  endif()
endif()

if(FTDQMC_ENABLE_SLEPC)
  # Try SLEPc CMake config first. If missing, fall back to pkg-config.
  find_package(SLEPc QUIET)
  if(NOT SLEPc_FOUND)
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
      if(DEFINED ENV{SLEPC_DIR} AND DEFINED ENV{PETSC_ARCH})
        set(ENV{PKG_CONFIG_PATH} "$ENV{SLEPC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
      endif()
      if(DEFINED ENV{PETSC_DIR} AND DEFINED ENV{PETSC_ARCH})
        set(ENV{PKG_CONFIG_PATH} "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
      endif()
      pkg_check_modules(SLEPC QUIET SLEPc slepc)
      if(SLEPC_FOUND)
        add_library(SLEPc::SLEPc INTERFACE IMPORTED)
        target_include_directories(SLEPc::SLEPc INTERFACE ${SLEPC_INCLUDE_DIRS})
        target_link_libraries(SLEPc::SLEPc INTERFACE ${SLEPC_LINK_LIBRARIES})
        set(SLEPc_FOUND TRUE)
      endif()
    endif()
  endif()
  if(NOT SLEPc_FOUND AND NOT SLEPC_FOUND)
    message(FATAL_ERROR "SLEPc not found. Set SLEPc_DIR or PETSC_DIR/PETSC_ARCH for pkg-config fallback.")
  endif()
endif()

# Build metadata
execute_process(
  COMMAND git rev-parse --short HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE FTDQMC_GIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)
if(NOT FTDQMC_GIT_HASH)
  set(FTDQMC_GIT_HASH "unknown")
endif()

string(TIMESTAMP FTDQMC_BUILD_TIME "%Y-%m-%d %H:%M:%S")

set(FTDQMC_COMPILER_ID "${CMAKE_CXX_COMPILER_ID}")
set(FTDQMC_COMPILER_VERSION "${CMAKE_CXX_COMPILER_VERSION}")

message(STATUS "FTDQMC-Studio build configuration:")
message(STATUS "  C++ compiler       : ${FTDQMC_COMPILER_ID} ${FTDQMC_COMPILER_VERSION}")
message(STATUS "  MPI                : ${FTDQMC_ENABLE_MPI}")
message(STATUS "  PETSc              : ${FTDQMC_ENABLE_PETSC}")
message(STATUS "  SLEPc              : ${FTDQMC_ENABLE_SLEPC}")
message(STATUS "  HDF5               : ${FTDQMC_ENABLE_HDF5}")
message(STATUS "  LINALG backend     : ${FTDQMC_LINALG_BACKEND}")

add_library(ftdqmc_core
  src/core/Logging.cpp
  src/config/Config.cpp
)

target_include_directories(ftdqmc_core PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(ftdqmc_core PUBLIC tomlplusplus::tomlplusplus)

if(FTDQMC_ENABLE_MPI)
  target_link_libraries(ftdqmc_core PUBLIC MPI::MPI_CXX)
  target_compile_definitions(ftdqmc_core PUBLIC FTDQMC_ENABLE_MPI=1)
else()
  target_compile_definitions(ftdqmc_core PUBLIC FTDQMC_ENABLE_MPI=0)
endif()

if(FTDQMC_ENABLE_PETSC)
  target_link_libraries(ftdqmc_core PUBLIC PETSc::PETSc)
  target_compile_definitions(ftdqmc_core PUBLIC FTDQMC_ENABLE_PETSC=1)
else()
  target_compile_definitions(ftdqmc_core PUBLIC FTDQMC_ENABLE_PETSC=0)
endif()

if(FTDQMC_ENABLE_SLEPC)
  target_link_libraries(ftdqmc_core PUBLIC SLEPc::SLEPc)
  target_compile_definitions(ftdqmc_core PUBLIC FTDQMC_ENABLE_SLEPC=1)
else()
  target_compile_definitions(ftdqmc_core PUBLIC FTDQMC_ENABLE_SLEPC=0)
endif()

if(FTDQMC_ENABLE_HDF5)
  target_compile_definitions(ftdqmc_core PUBLIC FTDQMC_ENABLE_HDF5=1)
else()
  target_compile_definitions(ftdqmc_core PUBLIC FTDQMC_ENABLE_HDF5=0)
endif()

target_compile_definitions(ftdqmc_core PUBLIC
  FTDQMC_GIT_HASH="${FTDQMC_GIT_HASH}"
  FTDQMC_BUILD_TIME="${FTDQMC_BUILD_TIME}"
  FTDQMC_COMPILER_ID="${FTDQMC_COMPILER_ID}"
  FTDQMC_COMPILER_VERSION="${FTDQMC_COMPILER_VERSION}"
  FTDQMC_LINALG_BACKEND="${FTDQMC_LINALG_BACKEND}"
)

add_executable(ftdqmc_studio apps/ftdqmc_studio/main.cpp)

target_link_libraries(ftdqmc_studio PRIVATE ftdqmc_core)

if(FTDQMC_ENABLE_MPI)
  target_compile_definitions(ftdqmc_studio PRIVATE FTDQMC_ENABLE_MPI=1)
else()
  target_compile_definitions(ftdqmc_studio PRIVATE FTDQMC_ENABLE_MPI=0)
endif()

if(FTDQMC_ENABLE_PETSC)
  target_compile_definitions(ftdqmc_studio PRIVATE FTDQMC_ENABLE_PETSC=1)
else()
  target_compile_definitions(ftdqmc_studio PRIVATE FTDQMC_ENABLE_PETSC=0)
endif()

if(FTDQMC_ENABLE_SLEPC)
  target_compile_definitions(ftdqmc_studio PRIVATE FTDQMC_ENABLE_SLEPC=1)
else()
  target_compile_definitions(ftdqmc_studio PRIVATE FTDQMC_ENABLE_SLEPC=0)
endif()

if(FTDQMC_ENABLE_HDF5)
  target_compile_definitions(ftdqmc_studio PRIVATE FTDQMC_ENABLE_HDF5=1)
else()
  target_compile_definitions(ftdqmc_studio PRIVATE FTDQMC_ENABLE_HDF5=0)
endif()
