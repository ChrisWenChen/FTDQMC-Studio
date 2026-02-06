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

LatticeType parse_lattice_type(const std::string& value) {
  const std::string v = to_lower(value);
  if (v == "honeycomb") {
    return LatticeType::Honeycomb;
  }
  if (v == "chain1d") {
    return LatticeType::Chain1D;
  }
  if (v == "rectangular") {
    return LatticeType::Rectangular;
  }
  if (v == "hexagonal") {
    return LatticeType::Hexagonal;
  }
  throw std::runtime_error("Invalid lattice.type: " + value);
}

BoundaryCondition parse_boundary_condition(const std::string& value) {
  const std::string v = to_lower(value);
  if (v == "obc") {
    return BoundaryCondition::OBC;
  }
  if (v == "pbc") {
    return BoundaryCondition::PBC;
  }
  if (v == "tbc") {
    return BoundaryCondition::TBC;
  }
  throw std::runtime_error("Invalid lattice.bc: " + value);
}

ModelType parse_model_type(const std::string& value) {
  const std::string v = to_lower(value);
  if (v == "hubbard") {
    return ModelType::Hubbard;
  }
  if (v == "haldane_hubbard") {
    return ModelType::HaldaneHubbard;
  }
  throw std::runtime_error("Invalid model.type: " + value);
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

  if (auto lattice = tbl["lattice"].as_table()) {
    if (auto type = (*lattice)["type"].value<std::string>()) {
      config.lattice.type = parse_lattice_type(*type);
    } else {
      throw std::runtime_error("Missing required [lattice].type");
    }
    if (auto Lx = (*lattice)["Lx"].value<int>()) {
      config.lattice.Lx = *Lx;
    }
    if (auto Ly = (*lattice)["Ly"].value<int>()) {
      config.lattice.Ly = *Ly;
    }
    if (auto bc = (*lattice)["bc"].value<std::string>()) {
      config.lattice.bc = parse_boundary_condition(*bc);
    }
    if (auto tx = (*lattice)["twist_theta_x"].value<double>()) {
      config.lattice.twist_theta_x = *tx;
    }
    if (auto ty = (*lattice)["twist_theta_y"].value<double>()) {
      config.lattice.twist_theta_y = *ty;
    }
  } else {
    throw std::runtime_error("Missing [lattice] table");
  }

  if (auto model = tbl["model"].as_table()) {
    if (auto type = (*model)["type"].value<std::string>()) {
      config.model.type = parse_model_type(*type);
    } else {
      throw std::runtime_error("Missing required [model].type");
    }
    if (auto U = (*model)["U"].value<double>()) {
      config.model.U = *U;
    }
    if (auto mu = (*model)["mu"].value<double>()) {
      config.model.mu = *mu;
    }
    if (auto t = (*model)["t"].value<double>()) {
      config.model.t = *t;
    }
    if (auto t2 = (*model)["t2"].value<double>()) {
      config.model.t2 = *t2;
    }
    if (auto phi = (*model)["phi"].value<double>()) {
      config.model.phi = *phi;
    }
  } else {
    throw std::runtime_error("Missing [model] table");
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
  if (config.lattice.Lx <= 0 || config.lattice.Ly <= 0) {
    throw std::runtime_error("[lattice].Lx and [lattice].Ly must be > 0");
  }
}

}  // namespace ftdqmc::config
