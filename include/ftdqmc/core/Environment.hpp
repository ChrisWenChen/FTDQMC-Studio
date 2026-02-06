#pragma once

#include <string>

namespace ftdqmc::core {

struct Environment {
  int rank = 0;
  int size = 1;
  std::string input_path;
  std::string out_dir;
};

}  // namespace ftdqmc::core
