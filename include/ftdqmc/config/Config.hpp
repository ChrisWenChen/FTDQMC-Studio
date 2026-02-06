#pragma once

#include <string>
#include <vector>

namespace ftdqmc::config {

enum class WorkflowMode { Spectrum, MC };

enum class LatticeType { Honeycomb, Chain1D, Rectangular, Hexagonal };

enum class BoundaryCondition { OBC, PBC, TBC };

enum class ModelType { Hubbard, HaldaneHubbard };

enum class TaskKind { Eigs, Svd, Measure, Unknown };

enum class TargetKind { M_mat, Bchain, Bprod, Unknown };

struct LatticeSpec {
  LatticeType type = LatticeType::Honeycomb;
  int Lx = 0;
  int Ly = 0;
  BoundaryCondition bc = BoundaryCondition::PBC;
  double twist_theta_x = 0.0;
  double twist_theta_y = 0.0;
};

struct ModelSpec {
  ModelType type = ModelType::Hubbard;
  double U = 0.0;
  double mu = 0.0;
  double t1 = 1.0;
  double t2 = 0.0;
  double phi = 0.0;
  double staggered_mass = 0.0;
};

struct TaskSpec {
  TaskKind kind = TaskKind::Unknown;
  TargetKind target = TargetKind::Unknown;
  std::string backend;
  int nev = 0;
};

struct Config {
  WorkflowMode workflow_mode = WorkflowMode::Spectrum;
  std::string run_name = "default_run";
  LatticeSpec lattice;
  ModelSpec model;
  std::string out_dir = "./out";
  bool per_rank_logs = true;
  std::vector<TaskSpec> tasks;
};

Config load_config_from_toml(const std::string& path);
void validate_config(const Config& config);

}  // namespace ftdqmc::config
