#pragma once

#include <fstream>
#include <string>

#include "ftdqmc/core/Environment.hpp"

namespace ftdqmc::core {

class Logger {
 public:
  Logger() = default;

  void open(const Environment& env, bool per_rank_logs);
  void info(const std::string& message);

 private:
  bool to_stdout_ = false;
  std::ofstream rank_file_;
  std::ofstream overview_file_;
};

Logger& log();

}  // namespace ftdqmc::core
