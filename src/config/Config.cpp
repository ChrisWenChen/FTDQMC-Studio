#include "ftdqmc/config/Config.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>

#include <toml++/toml.h>

namespace ftdqmc::config {

namespace {
std::string to_lower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return value;
}

WorkflowMode parse_workflow_mode(const std::string& value) {
  const std::string v = to_lower(value);
  if (v == "spectrum") {
    return WorkflowMode::Spectrum;
  }
  if (v == "mc") {
    return WorkflowMode::MC;
  }
  throw std::runtime_error("Invalid workflow.mode: " + value);
}

TaskKind parse_task_kind(const std::string& value) {
  const std::string v = to_lower(value);
  if (v == "eigs") {
    return TaskKind::Eigs;
  }
  if (v == "svd") {
    return TaskKind::Svd;
  }
  if (v == "measure") {
    return TaskKind::Measure;
  }
  return TaskKind::Unknown;
}

TargetKind parse_target_kind(const std::string& value) {
  const std::string v = to_lower(value);
  if (v == "m_mat") {
    return TargetKind::M_mat;
  }
  if (v == "bchain") {
    return TargetKind::Bchain;
  }
  if (v == "bprod") {
    return TargetKind::Bprod;
  }
  return TargetKind::Unknown;
}

}  // namespace

Config load_config_from_toml(const std::string& path) {
  Config config;

  toml::table tbl = toml::parse_file(path);

  if (auto workflow = tbl["workflow"].as_table()) {
    if (auto mode = (*workflow)["mode"].value<std::string>()) {
      config.workflow_mode = parse_workflow_mode(*mode);
    } else {
      throw std::runtime_error("Missing required [workflow].mode");
    }
    if (auto run_name = (*workflow)["run_name"].value<std::string>()) {
      config.run_name = *run_name;
    }
  } else {
    throw std::runtime_error("Missing [workflow] table");
  }

  if (auto io = tbl["io"].as_table()) {
    if (auto out_dir = (*io)["out_dir"].value<std::string>()) {
      config.out_dir = *out_dir;
    }
    if (auto per_rank_logs = (*io)["per_rank_logs"].value<bool>()) {
      config.per_rank_logs = *per_rank_logs;
    }
  }

  if (auto tasks = tbl["tasks"].as_array()) {
    for (auto&& node : *tasks) {
      auto task_table = node.as_table();
      if (!task_table) {
        continue;
      }
      TaskSpec task;
      if (auto kind = (*task_table)["kind"].value<std::string>()) {
        task.kind = parse_task_kind(*kind);
      }
      if (auto target = (*task_table)["target"].value<std::string>()) {
        task.target = parse_target_kind(*target);
      }
      if (auto backend = (*task_table)["backend"].value<std::string>()) {
        task.backend = *backend;
      }
      if (auto nev = (*task_table)["nev"].value<int>()) {
        task.nev = *nev;
      }
      config.tasks.push_back(task);
    }
  }

  return config;
}

void validate_config(const Config& config) {
  if (config.out_dir.empty()) {
    throw std::runtime_error("[io].out_dir must not be empty");
  }
  (void)config;
}

}  // namespace ftdqmc::config
