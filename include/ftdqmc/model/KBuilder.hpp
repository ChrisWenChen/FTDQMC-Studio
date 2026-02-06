#pragma once

#include "ftdqmc/linalg/DenseMatrix.hpp"
#include "ftdqmc/lattice/Honeycomb.hpp"
#include "ftdqmc/model/Model.hpp"

namespace ftdqmc::model {

linalg::DenseMatrix build_dense(const lattice::HoneycombLattice& lattice, const ModelSpec& model);

}  // namespace ftdqmc::model
