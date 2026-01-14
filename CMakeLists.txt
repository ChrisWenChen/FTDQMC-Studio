cmake_minimum_required(VERSION 3.20)

project(FTDQMC_Studio LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# (Optional but useful) Print compiler info for debugging
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION})")

# Define a user-configurable cache variable to select the BLAS backend with a fixed set of options.
set(FTDQMC_BLAS_BACKEND "auto" CACHE STRING "BLAS backend: auto|mkl|openblas|accelerate")
set_property(CACHE FTDQMC_BLAS_BACKEND PROPERTY STRINGS auto mkl openblas accelerate)

option(FTDQMC_ENABLE_MPI "Enable MPI" OFF)
option(FTDQMC_ENABLE_OPENMP "Enable OpenMP" OFF)
option(FTDQMC_ENABLE_PETSC "Enable PETSc" OFF)
option(FTDQMC_ENABLE_SLEPC "Enable SLEPc (requires PETSc)" OFF)

# Automatically select the best BLAS library for the current Operating System
if(FTDQMC_BLAS_BACKEND STREQUAL "auto")
  if(APPLE)
    set(FTDQMC_BLAS_BACKEND "accelerate") # Choose Accelerate framework for macOS
  else()
    set(FTDQMC_BLAS_BACKEND "mkl")        # Choose Intel MKL for Linux systems
  endif()
endif()

# Define the executable program and its source files
add_executable(ftdqmc_demo src/main.cpp)

# ---------- MPI Configuration ----------
if(FTDQMC_ENABLE_MPI)
  # REQUIRED makes CMake fail fast with a standard diagnostic message if MPI is missing
  find_package(MPI REQUIRED)
  target_link_libraries(ftdqmc_demo PRIVATE MPI::MPI_CXX)
  target_compile_definitions(ftdqmc_demo PRIVATE FTDQMC_USE_MPI=1)
endif()

# ---------- OpenMP Configuration ----------
if(FTDQMC_ENABLE_OPENMP)
  # REQUIRED makes CMake fail fast with a standard diagnostic message if OpenMP is missing
  find_package(OpenMP REQUIRED)
  target_link_libraries(ftdqmc_demo PRIVATE OpenMP::OpenMP_CXX)
  target_compile_definitions(ftdqmc_demo PRIVATE FTDQMC_USE_OPENMP=1)
endif()

# ---------- BLAS/LAPACK Configuration ----------
# NOTE: We link both BLAS and LAPACK because future kernels will need LAPACK (zgesvd/zheev/etc.)
# Use mutually-exclusive selection via elseif for safety and clarity.
if(FTDQMC_BLAS_BACKEND STREQUAL "accelerate")
  # ---------- Accelerate Configuration (macOS) ----------
  if(NOT APPLE)
    message(FATAL_ERROR "FTDQMC Error: Accelerate framework is only available on macOS platforms.") # Ensure the backend matches the operating system
  endif()

  # REQUIRED gives a standard CMake error if not found
  find_library(ACCELERATE_FRAMEWORK Accelerate REQUIRED)
  target_link_libraries(ftdqmc_demo PRIVATE ${ACCELERATE_FRAMEWORK})
  target_compile_definitions(ftdqmc_demo PRIVATE FTDQMC_BLAS_ACCELERATE=1)

elseif(FTDQMC_BLAS_BACKEND STREQUAL "mkl")
  # ---------- MKL Configuration (Intel) ----------
  # Prefer MKL CONFIG packages if available (common with oneMKL installs).
  find_package(MKL CONFIG QUIET)

  if(MKL_FOUND)
    target_link_libraries(ftdqmc_demo PRIVATE MKL::MKL)
    target_compile_definitions(ftdqmc_demo PRIVATE FTDQMC_BLAS_MKL=1)
  else()
    # Fallback: try FindBLAS/FindLAPACK with Intel vendor.
    # This typically works if environment is sourced (e.g., oneAPI setvars.sh) or MKL is installed system-wide.
    set(BLA_VENDOR Intel10_64lp)
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    target_link_libraries(ftdqmc_demo PRIVATE BLAS::BLAS LAPACK::LAPACK)
    target_compile_definitions(ftdqmc_demo PRIVATE FTDQMC_BLAS_MKL=1)
  endif()

elseif(FTDQMC_BLAS_BACKEND STREQUAL "openblas")
  # ---------- OpenBLAS Configuration ----------
  set(BLA_VENDOR OpenBLAS)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  target_link_libraries(ftdqmc_demo PRIVATE BLAS::BLAS LAPACK::LAPACK)
  target_compile_definitions(ftdqmc_demo PRIVATE FTDQMC_BLAS_OPENBLAS=1)

else()
  message(FATAL_ERROR "FTDQMC Error: Unknown FTDQMC_BLAS_BACKEND='${FTDQMC_BLAS_BACKEND}'. Use auto|mkl|openblas|accelerate.")
endif()

# ---------- PETSc Configuration ----------
if(FTDQMC_ENABLE_PETSC)
  # REQUIRED gives a standard CMake error if pkg-config is missing
  find_package(PkgConfig REQUIRED)

  # REQUIRED makes CMake stop with a standard error if PETSc is not found via pkg-config.
  # You can control the search path by setting PKG_CONFIG_PATH (often pointing to .../lib/pkgconfig).
  pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET petsc)

  target_link_libraries(ftdqmc_demo PRIVATE PkgConfig::PETSC)
  target_compile_definitions(ftdqmc_demo PRIVATE FTDQMC_USE_PETSC=1)

  # ---------- SLEPc Configuration ----------
  if(FTDQMC_ENABLE_SLEPC)
    # SLEPc requires PETSc; REQUIRED ensures clean failure if missing.
    pkg_check_modules(SLEPC REQUIRED IMPORTED_TARGET slepc)

    target_link_libraries(ftdqmc_demo PRIVATE PkgConfig::SLEPC)
    target_compile_definitions(ftdqmc_demo PRIVATE FTDQMC_USE_SLEPC=1)
  endif()
endif()

# ---------- Testing Configuration ----------
# include(CTest) sets up standard test variables (e.g., BUILD_TESTING) and also calls enable_testing().
include(CTest)
enable_testing() # Core command to enable the CTest framework

option(FTDQMC_ENABLE_GTEST "Enable GoogleTest unit tests" ON)

if(FTDQMC_ENABLE_GTEST AND BUILD_TESTING)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  )
  # Ensures the googletest targets (GTest::gtest_main, etc.) are actually created.
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  add_executable(test_smoke tests/test_smoke.cpp)
  target_link_libraries(test_smoke PRIVATE GTest::gtest_main)

  # Validate the test target exists (TARGET is the correct way to check, not "if(NOT test_smoke)")
  if(NOT TARGET test_smoke)
    message(FATAL_ERROR "FTDQMC Error: Failed to create test executable!") # Check for test target creation failure
  endif()

  include(GoogleTest)
  gtest_discover_tests(test_smoke) # Automatically find and register GTest units to CTest
endif()
