#pragma once

#include <string>

namespace ftdqmc::model {

enum class ModelType { Hubbard, HaldaneHubbard };

struct HubbardParams {
  double U = 0.0;
  double mu = 0.0;
  double t = 1.0;
};

struct HaldaneParams {
  double U = 0.0;
  double mu = 0.0;
  double t = 1.0;
  double t2 = 0.0;
  double phi = 0.0;
};

struct ModelSpec {
  ModelType type = ModelType::Hubbard;
  HubbardParams hubbard;
  HaldaneParams haldane;
};

}  // namespace ftdqmc::model
