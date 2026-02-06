#include "ftdqmc/model/KBuilder.hpp"

#include <complex>

namespace ftdqmc::model {

namespace {
void add_hopping(linalg::DenseMatrix& K, int i, int j, std::complex<double> t) {
  K(i, j) += t;
  K(j, i) += std::conj(t);
}
}

linalg::DenseMatrix build_dense(const lattice::HoneycombLattice& lattice, const ModelSpec& model) {
  const int Ns = lattice.Ns();
  linalg::DenseMatrix K = linalg::DenseMatrix::zeros(Ns, Ns);

  if (model.type == ModelType::Hubbard) {
    const auto& p = model.hubbard;
    for (const auto& bond : lattice.nn_bonds()) {
      if (bond.kind != lattice::BondKind::NN) {
        continue;
      }
      add_hopping(K, bond.i, bond.j, -p.t1);
    }
    for (int site = 0; site < Ns; ++site) {
      K(site, site) += -p.mu;
    }
    return K;
  }

  const auto& p = model.haldane;
  for (const auto& bond : lattice.nn_bonds()) {
    if (bond.kind != lattice::BondKind::NN) {
      continue;
    }
    add_hopping(K, bond.i, bond.j, -p.t1);
  }

  for (const auto& bond : lattice.nnn_bonds()) {
    if (bond.kind != lattice::BondKind::NNN) {
      continue;
    }
    const std::complex<double> phase = std::exp(std::complex<double>(0.0, bond.orientation * p.phi));
    const std::complex<double> t = -p.t2 * phase;
    add_hopping(K, bond.i, bond.j, t);
  }

  for (int site = 0; site < Ns; ++site) {
    const auto sub = lattice.sublattice(site);
    const double m = (sub == lattice::Sublattice::A) ? p.staggered_mass : -p.staggered_mass;
    K(site, site) += -p.mu + m;
  }

  return K;
}

}  // namespace ftdqmc::model
