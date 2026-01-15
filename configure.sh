#!/usr/bin/env bash
set -euo pipefail

# ===============================
# Defaults (PETSc-style philosophy)
# ===============================

OS="$(uname)"
OS_LOWER="$(echo "$OS" | tr '[:upper:]' '[:lower:]')"
BUILD_TYPE="Release"

# Feature toggles
MPI="OFF"
OPENMP="OFF"
PETSC="OFF"
SLEPC="OFF"
GTEST="ON"

# BLAS backend
BLAS="auto"

# Compiler defaults by OS
if [[ "$OS" == "Darwin" ]]; then
  COMPILER="clang++"
  BLAS="accelerate"
else
  COMPILER="icpx"
  BLAS="mkl"
fi

# PETSc/SLEPc paths (optional)
PETSC_DIR_ARG=""
SLEPC_DIR_ARG=""
PETSC_ARCH_ARG="${PETSC_ARCH:-}"

# Build directory
BUILD_DIR="build/${OS_LOWER}"

# ===============================
# Help
# ===============================
usage() {
  cat <<EOF
Usage: ./configure.sh [options]

Compiler / backend:
  --compiler=<g++|clang++|icpx|mpic++|mpiicpx>
  --blas=<auto|mkl|openblas|accelerate>

Features:
  --mpi=on|off
  --openmp=on|off
  --petsc=on|off
  --slepc=on|off        (requires PETSc)
  --gtest=on|off

PETSc / SLEPc paths (optional):
  --petsc-dir=<path>
  --slepc-dir=<path>
  --petsc-arch=<arch>

Build:
  --type=<Release|Debug>
  --build-dir=<path>

Examples:
  ./configure.sh
  ./configure.sh --compiler=g++ --blas=openblas
  ./configure.sh --mpi=on --petsc=on --slepc=on --petsc-dir=\$PETSC_DIR --petsc-arch=\$PETSC_ARCH --slepc-dir=\$SLEPC_DIR
EOF
}

# ===============================
# Parse arguments
# ===============================
for arg in "$@"; do
  case "$arg" in
    --compiler=*)   COMPILER="${arg#*=}" ;;
    --blas=*)       BLAS="${arg#*=}" ;;
    --mpi=on)       MPI="ON" ;;
    --mpi=off)      MPI="OFF" ;;
    --openmp=on)    OPENMP="ON" ;;
    --openmp=off)   OPENMP="OFF" ;;
    --petsc=on)     PETSC="ON" ;;
    --petsc=off)    PETSC="OFF" ;;
    --slepc=on)     SLEPC="ON" ;;
    --slepc=off)    SLEPC="OFF" ;;
    --gtest=on)     GTEST="ON" ;;
    --gtest=off)    GTEST="OFF" ;;
    --petsc-dir=*)  PETSC_DIR_ARG="${arg#*=}" ;;
    --slepc-dir=*)  SLEPC_DIR_ARG="${arg#*=}" ;;
    --petsc-arch=*) PETSC_ARCH_ARG="${arg#*=}" ;;
    --type=*)       BUILD_TYPE="${arg#*=}" ;;
    --build-dir=*)  BUILD_DIR="${arg#*=}" ;;
    -h|--help)      usage; exit 0 ;;
    *)
      echo "Unknown option: $arg"
      usage
      exit 1
      ;;
  esac
done

# ===============================
# Sanity checks & fallbacks
# ===============================

# Compiler fallback
if ! command -v "$COMPILER" >/dev/null 2>&1; then
  echo "[configure] Compiler '$COMPILER' not found -> falling back to g++"
  COMPILER="g++"
fi

# accelerate only on macOS
if [[ "$OS" != "Darwin" && "$BLAS" == "accelerate" ]]; then
  echo "[configure] Error: Accelerate is only available on macOS"
  exit 1
fi

# slepc implies petsc
if [[ "$SLEPC" == "ON" && "$PETSC" == "OFF" ]]; then
  echo "[configure] slepc=on requires petsc -> enabling PETSc"
  PETSC="ON"
fi

# ===============================
# Handle PETSc / SLEPc pkg-config paths
# ===============================

PKG_PATH_EXTRA=""

if [[ -n "$PETSC_DIR_ARG" ]]; then
  if [[ -n "$PETSC_ARCH_ARG" && -d "$PETSC_DIR_ARG/$PETSC_ARCH_ARG/lib/pkgconfig" ]]; then
    PKG_PATH_EXTRA="$PETSC_DIR_ARG/$PETSC_ARCH_ARG/lib/pkgconfig"
  elif [[ -d "$PETSC_DIR_ARG/lib/pkgconfig" ]]; then
    PKG_PATH_EXTRA="$PETSC_DIR_ARG/lib/pkgconfig"
  fi
fi

if [[ -n "$SLEPC_DIR_ARG" ]]; then
  if [[ -n "$PETSC_ARCH_ARG" && -d "$SLEPC_DIR_ARG/$PETSC_ARCH_ARG/lib/pkgconfig" ]]; then
    PKG_PATH_EXTRA="$PKG_PATH_EXTRA:$SLEPC_DIR_ARG/$PETSC_ARCH_ARG/lib/pkgconfig"
  elif [[ -d "$SLEPC_DIR_ARG/lib/pkgconfig" ]]; then
    PKG_PATH_EXTRA="$PKG_PATH_EXTRA:$SLEPC_DIR_ARG/lib/pkgconfig"
  fi
fi

if [[ -n "$PKG_PATH_EXTRA" ]]; then
  export PKG_CONFIG_PATH="$PKG_PATH_EXTRA:${PKG_CONFIG_PATH:-}"
  echo "[configure] PKG_CONFIG_PATH += $PKG_PATH_EXTRA"
fi

# ===============================
# Summary
# ===============================
echo "----------------------------------------"
echo "FTDQMC-Studio configure summary"
echo "  OS            = $OS"
echo "  Build type    = $BUILD_TYPE"
echo "  CXX compiler  = $COMPILER"
echo "  BLAS backend  = $BLAS"
echo "  MPI           = $MPI"
echo "  OpenMP        = $OPENMP"
echo "  PETSc         = $PETSC"
echo "  SLEPc         = $SLEPC"
echo "  GTest         = $GTEST"
echo "  Build dir     = $BUILD_DIR"
echo "----------------------------------------"

# ===============================
# Run CMake configure
# ===============================
cmake -S . -B "$BUILD_DIR" \
  -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
  -DCMAKE_CXX_COMPILER="$COMPILER" \
  -DFTDQMC_BLAS_BACKEND="$BLAS" \
  -DFTDQMC_ENABLE_MPI="$MPI" \
  -DFTDQMC_ENABLE_OPENMP="$OPENMP" \
  -DFTDQMC_ENABLE_PETSC="$PETSC" \
  -DFTDQMC_ENABLE_SLEPC="$SLEPC" \
  -DFTDQMC_ENABLE_GTEST="$GTEST"

echo
echo "[configure] Configuration completed successfully."
echo "Next steps:"
echo "  cmake --build $BUILD_DIR -j"
echo "  $BUILD_DIR/ftdqmc_demo"
echo "  ctest --test-dir $BUILD_DIR --output-on-failure"
