#include "ftdqmc/lattice/chain1d.hpp"

namespace ftdqmc {

Chain1D::Chain1D(int Lx, bool pbc_x)
  : Lx_(Lx), pbc_x_(pbc_x) {

  sites_.reserve(Lx_);
  for (int x = 0; x < Lx_; ++x) {
    sites_.push_back(Site{x, 0, 0});
  }

  // NN: i -> i+1 (directed)
  for (int x = 0; x < Lx_; ++x) {
    int i = idx(x);
    int nx = x + 1;
    int dx = 0;
    if (nx < Lx_ || pbc_x_) {
      if (nx >= Lx_) { 
        nx = wrap(nx, Lx_);
        dx = +1; 
      }
      nn_.push_back(Bond{i, idx(nx), dx, 0});
    }
  }
}

} // namespace ftdqmc
