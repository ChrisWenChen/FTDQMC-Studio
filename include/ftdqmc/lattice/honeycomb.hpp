#pragma once
#include "ftdqmc/lattice/lattice.hpp"

namespace ftdqmc {

class HoneycombLattice final : public Lattice {
public:
  HoneycombLattice(int Lx, int Ly, bool pbc_x=true, bool pbc_y=true, bool include_nnn=false);

private:
  int Lx_, Ly_;
  bool pbc_x_, pbc_y_;

  // index: (x,y,sub) -> [0, 2*Lx*Ly)
  int idx(int x, int y, int sub) const { return sub + 2 * (x + Lx_ * y); }
  int wrap(int a, int L) const { return (a % L + L) % L; }
};

} // namespace ftdqmc
