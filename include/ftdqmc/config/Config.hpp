#pragma once

#include <string>
#include <vector>

namespace ftdqmc::config {

enum class WorkflowMode { Spectrum, MC };

enum class TaskKind { Eigs, Svd, Measure, Unknown };

enum class TargetKind { M_mat, Bchain, Bprod, Unknown };

struct TaskSpec {
  TaskKind kind = TaskKind::Unknown;
  TargetKind target = TargetKind::Unknown;
  std::string backend;
  int nev = 0;
};

struct Config {
  WorkflowMode workflow_mode = WorkflowMode::Spectrum;
  std::string run_name = "default_run";
  std::string out_dir = "./out";
  bool per_rank_logs = true;
  std::vector<TaskSpec> tasks;
};

Config load_config_from_toml(const std::string& path);
void validate_config(const Config& config);

}  // namespace ftdqmc::config
