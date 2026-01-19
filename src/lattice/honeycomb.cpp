#include "ftdqmc/lattice/honeycomb.hpp"

namespace ftdqmc {

HoneycombLattice::HoneycombLattice(int Lx, int Ly, bool pbc_x, bool pbc_y, bool include_nnn)
  : Lx_(Lx), Ly_(Ly), pbc_x_(pbc_x), pbc_y_(pbc_y) {

  // Initialize sites
  sites_.reserve(2 * Lx_ * Ly_);
  for (int y = 0; y < Ly_; ++y) {
    for (int x = 0; x < Lx_; ++x) {
      sites_.push_back(Site{x, y, 0}); // A
      sites_.push_back(Site{x, y, 1}); // B
    }
  }

  // NN bonds: A(x,y) -> B neighbors (3) in brick-wall convention:
  //   B(x,y), B(x-1,y), B(x,y-1)
  for (int y = 0; y < Ly_; ++y) {
    for (int x = 0; x < Lx_; ++x) {
      int a = idx(x, y, 0);

      // Inside the unit cell
      nn_.push_back(Bond{a, idx(x, y, 1), 0, 0}); // B(x,y)

      // B(x-1,y)
      {
        int nx = x - 1;
        int dx = 0;
        if (nx >= 0 || pbc_x_) { 
          if (nx < 0) {
            nx = wrap(nx, Lx_);
            dx = -1;
          }
          nn_.push_back(Bond{a, idx(nx, y, 1), dx, 0});
        }
      }

      // B(x,y-1)
      {
        int ny = y - 1;
        int dy = 0;
        if (ny >= 0 || pbc_y_) { 
          if (ny < 0) {
            ny = wrap(ny, Ly_);
            dy = -1;
          }
          nn_.push_back(Bond{a, idx(x, ny, 1), 0, dy});
        }
      }
    }
  }

  if (include_nnn) {
    // NNN bonds: same sublattice (A-A and B-B), 6 per site (store directed)
    const int dirs[6][2] = {
      {+1, 0}, {0, +1}, {+1, +1},
      {-1, 0}, {0, -1}, {-1, -1}
    };

    for (int y = 0; y < Ly_; ++y) {
      for (int x = 0; x < Lx_; ++x) {
        for (int sub = 0; sub < 2; ++sub) {
          int i = idx(x, y, sub);
          for (auto &d : dirs) {
            int nx = x + d[0], ny = y + d[1];
            int dx = 0, dy = 0;

            bool x_valid = (nx >= 0 && nx < Lx_) || pbc_x_;
            bool y_valid = (ny >= 0 && ny < Ly_) || pbc_y_;

            if (x_valid && y_valid) {
              if (nx < 0 || nx >= Lx_) {
                dx = (nx < 0) ? -1 : +1;
                nx = wrap(nx, Lx_);
              }
              if (ny < 0 || ny >= Ly_) {
                dy = (ny < 0) ? -1 : +1;
                ny = wrap(ny, Ly_);
              }
              nnn_.push_back(Bond{i, idx(nx, ny, sub), dx, dy});
            }
          }
        }
      }
    }
  }
}

} // namespace ftdqmc
