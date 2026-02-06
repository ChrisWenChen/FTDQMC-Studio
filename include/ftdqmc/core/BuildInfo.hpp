#pragma once

#include <string>

namespace ftdqmc::core {

inline std::string git_hash() {
  return std::string(FTDQMC_GIT_HASH);
}

inline std::string build_time() {
  return std::string(FTDQMC_BUILD_TIME);
}

inline std::string compiler_info() {
  return std::string(FTDQMC_COMPILER_ID) + " " + std::string(FTDQMC_COMPILER_VERSION);
}

inline std::string linalg_backend() {
  return std::string(FTDQMC_LINALG_BACKEND);
}

inline std::string backend_flags() {
  std::string flags;
  flags += "MPI=" + std::string(FTDQMC_ENABLE_MPI ? "ON" : "OFF");
  flags += ", PETSc=" + std::string(FTDQMC_ENABLE_PETSC ? "ON" : "OFF");
  flags += ", SLEPc=" + std::string(FTDQMC_ENABLE_SLEPC ? "ON" : "OFF");
  flags += ", HDF5=" + std::string(FTDQMC_ENABLE_HDF5 ? "ON" : "OFF");
  flags += ", LINALG=" + linalg_backend();
  return flags;
}

}  // namespace ftdqmc::core
