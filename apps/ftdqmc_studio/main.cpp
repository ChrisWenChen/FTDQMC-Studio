#include <filesystem>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <optional>
#include <string>

#include "ftdqmc/Config.hpp"
#include "ftdqmc/core/BuildInfo.hpp"
#include "ftdqmc/core/Environment.hpp"
#include "ftdqmc/core/Logging.hpp"
#include "ftdqmc/linalg/DenseMatrix.hpp"
#include "ftdqmc/lattice/Honeycomb.hpp"
#include "ftdqmc/model/KBuilder.hpp"
#include "ftdqmc/model/Model.hpp"

#if FTDQMC_ENABLE_PETSC
  #include <petsc.h>
#elif FTDQMC_ENABLE_MPI
  #include <mpi.h>
#endif

namespace {
const char* kUsage = "Usage: ftdqmc_studio <input_file> [PETSc/SLEPc options...]";

std::optional<int> read_env_int(const char* name) {
  const char* value = std::getenv(name);
  if (!value || value[0] == '\0') {
    return std::nullopt;
  }
  try {
    return std::stoi(value);
  } catch (...) {
    return std::nullopt;
  }
}

bool running_under_mpi(size_t* out_size) {
  const char* envs[] = {
      "OMPI_COMM_WORLD_SIZE",
      "PMI_SIZE",
      "PMIX_SIZE",
      "MV2_COMM_WORLD_SIZE",
      "MPI_LOCALNRANKS",
  };
  for (const char* env : envs) {
    if (auto size = read_env_int(env)) {
      if (out_size) {
        *out_size = static_cast<size_t>(*size);
      }
      return true;
    }
  }
  return false;
}

std::string workflow_mode_to_string(ftdqmc::config::WorkflowMode mode) {
  switch (mode) {
    case ftdqmc::config::WorkflowMode::Spectrum:
      return "spectrum";
    case ftdqmc::config::WorkflowMode::MC:
      return "mc";
    default:
      return "unknown";
  }
}

std::string task_kind_to_string(ftdqmc::config::TaskKind kind) {
  switch (kind) {
    case ftdqmc::config::TaskKind::Eigs:
      return "eigs";
    case ftdqmc::config::TaskKind::Svd:
      return "svd";
    case ftdqmc::config::TaskKind::Measure:
      return "measure";
    default:
      return "unknown";
  }
}

std::string target_kind_to_string(ftdqmc::config::TargetKind target) {
  switch (target) {
    case ftdqmc::config::TargetKind::M_mat:
      return "M_mat";
    case ftdqmc::config::TargetKind::Bchain:
      return "Bchain";
    case ftdqmc::config::TargetKind::Bprod:
      return "Bprod";
    default:
      return "unknown";
  }
}

void log_task_summary(const ftdqmc::config::Config& config) {
  auto& logger = ftdqmc::core::log();
  logger.info("tasks.count=" + std::to_string(config.tasks.size()));
  for (size_t i = 0; i < config.tasks.size(); ++i) {
    const auto& task = config.tasks[i];
    std::string line = "tasks[" + std::to_string(i) + "]: kind=" +
                       task_kind_to_string(task.kind) + ", target=" +
                       target_kind_to_string(task.target);
    if (!task.backend.empty()) {
      line += ", backend=" + task.backend;
    }
    if (task.nev > 0) {
      line += ", nev=" + std::to_string(task.nev);
    }
    logger.info(line);
  }
}

ftdqmc::lattice::BoundaryCondition map_bc(ftdqmc::config::BoundaryCondition bc) {
  using BC = ftdqmc::config::BoundaryCondition;
  switch (bc) {
    case BC::OBC:
      return ftdqmc::lattice::BoundaryCondition::OBC;
    case BC::PBC:
      return ftdqmc::lattice::BoundaryCondition::PBC;
    case BC::TBC:
      return ftdqmc::lattice::BoundaryCondition::TBC;
    default:
      return ftdqmc::lattice::BoundaryCondition::PBC;
  }
}

ftdqmc::model::ModelSpec map_model(const ftdqmc::config::ModelSpec& in) {
  ftdqmc::model::ModelSpec out;
  if (in.type == ftdqmc::config::ModelType::HaldaneHubbard) {
    out.type = ftdqmc::model::ModelType::HaldaneHubbard;
    out.haldane.U = in.U;
    out.haldane.mu = in.mu;
    out.haldane.t1 = in.t1;
    out.haldane.t2 = in.t2;
    out.haldane.phi = in.phi;
    out.haldane.staggered_mass = in.staggered_mass;
  } else {
    out.type = ftdqmc::model::ModelType::Hubbard;
    out.hubbard.U = in.U;
    out.hubbard.mu = in.mu;
    out.hubbard.t1 = in.t1;
  }
  return out;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2 || argv[1][0] == '-') {
    std::cerr << kUsage << "\n";
    return 1;
  }

  const std::string input_path = argv[1];

#if !FTDQMC_ENABLE_MPI && !FTDQMC_ENABLE_PETSC
  size_t detected_size = 0;
  if (running_under_mpi(&detected_size) && detected_size > 1) {
    std::cerr << "Error: MPI is disabled at build time but mpirun was used (size="
              << detected_size << "). Rebuild with FTDQMC_ENABLE_MPI=ON.\n";
    return 2;
  }
#endif

#if FTDQMC_ENABLE_PETSC
  PetscInitialize(&argc, &argv, nullptr, nullptr);
#elif FTDQMC_ENABLE_MPI
  MPI_Init(&argc, &argv);
#endif

  int rank = 0;
  int size = 1;
#if FTDQMC_ENABLE_PETSC || FTDQMC_ENABLE_MPI
  int mpi_initialized = 0;
  MPI_Initialized(&mpi_initialized);
  if (mpi_initialized) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  }
#endif

  try {
    auto config = ftdqmc::config::load_config_from_toml(input_path);
    ftdqmc::config::validate_config(config);

    if (rank == 0) {
      std::filesystem::create_directories(config.out_dir);
    }
#if FTDQMC_ENABLE_PETSC || FTDQMC_ENABLE_MPI
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    if (mpi_initialized) {
      MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    ftdqmc::core::Environment env;
    env.rank = rank;
    env.size = size;
    env.input_path = input_path;
    env.out_dir = config.out_dir;

    auto& logger = ftdqmc::core::log();
    logger.open(env, config.per_rank_logs);

    logger.info("ftdqmc_studio start");
    logger.info("git_hash=" + ftdqmc::core::git_hash());
    logger.info("build_time=" + ftdqmc::core::build_time());
    logger.info("compiler=" + ftdqmc::core::compiler_info());
    logger.info("backend_flags=" + ftdqmc::core::backend_flags());
    logger.info("mpi_size=" + std::to_string(size) + ", mpi_rank=" + std::to_string(rank));
    logger.info("input_file=" + input_path);
    logger.info("out_dir=" + config.out_dir);
    logger.info("workflow.mode=" + workflow_mode_to_string(config.workflow_mode));
    logger.info("workflow.run_name=" + config.run_name);
    log_task_summary(config);

    if (config.lattice.type != ftdqmc::config::LatticeType::Honeycomb) {
      throw std::runtime_error("Milestone 1 supports only honeycomb lattice");
    }
    auto bc = map_bc(config.lattice.bc);
    ftdqmc::lattice::HoneycombLattice lattice(config.lattice.Lx, config.lattice.Ly, bc,
                                              config.lattice.twist_theta_x,
                                              config.lattice.twist_theta_y);
    if (config.lattice.bc == ftdqmc::config::BoundaryCondition::TBC) {
      logger.info("twist BC stub: theta_x=" + std::to_string(config.lattice.twist_theta_x) +
                  ", theta_y=" + std::to_string(config.lattice.twist_theta_y));
    }
    logger.info("lattice.Ns=" + std::to_string(lattice.Ns()) +
                ", Lx=" + std::to_string(lattice.Lx()) +
                ", Ly=" + std::to_string(lattice.Ly()));
    logger.info("lattice.nn_bonds=" + std::to_string(lattice.nn_bonds().size()) +
                ", lattice.nnn_bonds=" + std::to_string(lattice.nnn_bonds().size()));

    const auto model = map_model(config.model);
    ftdqmc::linalg::DenseMatrix K = ftdqmc::model::build_dense(lattice, model);
    const double herm_err = ftdqmc::linalg::hermitian_error_norm(K);
    const auto tr = ftdqmc::linalg::trace(K);
    logger.info("K.hermitian_error_norm=" + std::to_string(herm_err));
    logger.info("K.trace=" + std::to_string(tr.real()) + " + " + std::to_string(tr.imag()) + "i");

    logger.info("generate/read field: placeholder");
    logger.info("build M_mat^sigma: placeholder");
    logger.info("spectrum: placeholder");

    if (rank == 0) {
      std::filesystem::path status_path = std::filesystem::path(config.out_dir) / "status.txt";
      std::ofstream status_file(status_path, std::ios::out | std::ios::trunc);
      status_file << "status=ok\n";
      status_file << "note=milestone0_placeholder\n";

      std::filesystem::path summary_path = std::filesystem::path(config.out_dir) / "run_summary.txt";
      std::ofstream summary_file(summary_path, std::ios::out | std::ios::trunc);
      summary_file << "ftdqmc_studio run summary\n";
      summary_file << "input_file=" << input_path << "\n";
      summary_file << "out_dir=" << config.out_dir << "\n";
      summary_file << "workflow.mode=" << workflow_mode_to_string(config.workflow_mode) << "\n";
      summary_file << "tasks.count=" << config.tasks.size() << "\n";
      summary_file << "status=ok\n";
    }

    logger.info("write output: create dummy status file");
    logger.info("ftdqmc_studio done");
  } catch (const std::exception& ex) {
    if (rank == 0) {
      std::cerr << "Error: " << ex.what() << "\n";
    }
#if FTDQMC_ENABLE_PETSC
    PetscFinalize();
#elif FTDQMC_ENABLE_MPI
    MPI_Finalize();
#endif
    return 1;
  }

#if FTDQMC_ENABLE_PETSC
  PetscFinalize();
#elif FTDQMC_ENABLE_MPI
  MPI_Finalize();
#endif

  return 0;
}
