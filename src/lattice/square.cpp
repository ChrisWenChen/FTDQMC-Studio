#include "ftdqmc/lattice/square.hpp"

namespace ftdqmc {

SquareLattice::SquareLattice(int Lx, int Ly, bool pbc_x, bool pbc_y)
  : Lx_(Lx), Ly_(Ly), pbc_x_(pbc_x), pbc_y_(pbc_y) {

  sites_.reserve(Lx_ * Ly_);
  for (int y = 0; y < Ly_; ++y) {
    for (int x = 0; x < Lx_; ++x) {
      sites_.push_back(Site{x, y, 0});
    }
  }

  // NN: right and up (directed)
  for (int y = 0; y < Ly_; ++y) {
    for (int x = 0; x < Lx_; ++x) {
      int i = idx(x, y);

      // +x
      {
        int nx = x + 1;
        int dx = 0;
        if (nx < Lx_ || pbc_x_) {
            if (nx >= Lx_) {
                nx = wrap(nx, Lx_);
                dx = +1;
        }
        nn_.push_back(Bond{i, idx(nx, y), dx, 0});
        }
      }

      // +y
      {
        int ny = y + 1;
        int dy = 0;
        if (ny < Ly_ || pbc_y_) {
          if (ny >= Ly_) {
            ny = wrap(ny, Ly_);
            dy = +1;
        }
        nn_.push_back(Bond{i, idx(x, ny), 0, dy});
        }
      }
    }
  }
}

} // namespace ftdqmc
