#include "ftdqmc/core/Logging.hpp"

#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace ftdqmc::core {

namespace {
std::string timestamp_now() {
  const auto now = std::chrono::system_clock::now();
  const auto t = std::chrono::system_clock::to_time_t(now);
  std::tm tm{};
#if defined(_WIN32)
  localtime_s(&tm, &t);
#else
  localtime_r(&t, &tm);
#endif
  std::ostringstream oss;
  oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
  return oss.str();
}
}  // namespace

void Logger::open(const Environment& env, bool per_rank_logs) {
  to_stdout_ = (env.rank == 0);
  if (per_rank_logs) {
    std::filesystem::path rank_path = std::filesystem::path(env.out_dir) /
                                      ("log_rank" + std::to_string(env.rank) + ".txt");
    rank_file_.open(rank_path, std::ios::out | std::ios::app);
  }
  if (env.rank == 0) {
    std::filesystem::path overview_path = std::filesystem::path(env.out_dir) / "log_overview.txt";
    overview_file_.open(overview_path, std::ios::out | std::ios::app);
  }
}

void Logger::info(const std::string& message) {
  const std::string line = "[" + timestamp_now() + "] " + message + "\n";
  if (to_stdout_) {
    std::cout << line;
    std::cout.flush();
  }
  if (rank_file_.is_open()) {
    rank_file_ << line;
    rank_file_.flush();
  }
  if (overview_file_.is_open()) {
    overview_file_ << line;
    overview_file_.flush();
  }
}

Logger& log() {
  static Logger logger;
  return logger;
}

}  // namespace ftdqmc::core
