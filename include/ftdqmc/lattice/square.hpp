#pragma once
#include "ftdqmc/lattice/lattice.hpp"

namespace ftdqmc {

class SquareLattice final : public Lattice {
public:
  SquareLattice(int Lx, int Ly, bool pbc_x=true, bool pbc_y=true);

private:
  int Lx_, Ly_;
  bool pbc_x_, pbc_y_;

  int idx(int x, int y) const { return x + Lx_ * y; }
  int wrap(int a, int L) const { return (a % L + L) % L; }
};

} // namespace ftdqmc
