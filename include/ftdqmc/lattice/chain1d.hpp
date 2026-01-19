#pragma once
#include "ftdqmc/lattice/lattice.hpp"

namespace ftdqmc {

class Chain1D final : public Lattice {
public:
  Chain1D(int Lx, bool pbc_x=true);

private:
  int Lx_;
  bool pbc_x_;

  int idx(int x) const { return x; }
  int wrap(int a, int L) const { return (a % L + L) % L; }
};

} // namespace ftdqmc
