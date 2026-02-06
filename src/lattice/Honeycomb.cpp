#include "ftdqmc/lattice/Honeycomb.hpp"

#include <stdexcept>

namespace ftdqmc::lattice {

HoneycombLattice::HoneycombLattice(int Lx, int Ly, BoundaryCondition bc,
                                   double twist_theta_x, double twist_theta_y)
    : Lx_(Lx), Ly_(Ly), bc_(bc), twist_theta_x_(twist_theta_x), twist_theta_y_(twist_theta_y) {
  if (Lx_ <= 0 || Ly_ <= 0) {
    throw std::runtime_error("HoneycombLattice: Lx and Ly must be > 0");
  }
  ns_ = 2 * Lx_ * Ly_;
  build_bonds();
}

int HoneycombLattice::site_index(int x, int y, Sublattice sub) const {
  return 2 * (x + Lx_ * y) + static_cast<int>(sub);
}

bool HoneycombLattice::map_coords(int x, int y, int& xw, int& yw, int& dx_wrap, int& dy_wrap) const {
  dx_wrap = 0;
  dy_wrap = 0;
  if (bc_ == BoundaryCondition::OBC) {
    if (x < 0 || x >= Lx_ || y < 0 || y >= Ly_) {
      return false;
    }
    xw = x;
    yw = y;
    return true;
  }
  // PBC/TBC: wrap into [0, Lx-1], [0, Ly-1]
  xw = x % Lx_;
  yw = y % Ly_;
  if (xw < 0) xw += Lx_;
  if (yw < 0) yw += Ly_;
  dx_wrap = (x - xw) / Lx_;
  dy_wrap = (y - yw) / Ly_;
  return true;
}

Sublattice HoneycombLattice::sublattice(int site) const {
  int sub = site % 2;
  return sub == 0 ? Sublattice::A : Sublattice::B;
}

void HoneycombLattice::build_bonds() {
  nn_bonds_.clear();
  nnn_bonds_.clear();

  // NN bonds: for each A site, connect to three B neighbors
  const int nn_dx[3] = {0, -1, 0};
  const int nn_dy[3] = {0, 0, -1};

  for (int y = 0; y < Ly_; ++y) {
    for (int x = 0; x < Lx_; ++x) {
      const int iA = site_index(x, y, Sublattice::A);
      for (int k = 0; k < 3; ++k) {
        int x2 = x + nn_dx[k];
        int y2 = y + nn_dy[k];
        int xw = 0, yw = 0, dx_wrap = 0, dy_wrap = 0;
        if (!map_coords(x2, y2, xw, yw, dx_wrap, dy_wrap)) {
          continue;
        }
        const int jB = site_index(xw, yw, Sublattice::B);
        if (iA == jB) {
          continue;
        }
        Bond b;
        b.i = iA;
        b.j = jB;
        b.dx_wrap = dx_wrap;
        b.dy_wrap = dy_wrap;
        b.kind = BondKind::NN;
        b.orientation = 0;
        nn_bonds_.push_back(b);
      }
    }
  }

  // NNN bonds: for each site, connect to three same-sublattice neighbors
  // Displacements chosen to avoid double counting (each undirected pair once).
  const int nnn_dx[3] = {1, 0, -1};
  const int nnn_dy[3] = {0, -1, 1};
  // Orientation convention: sign depends on displacement direction.
  const int nnn_orient[3] = {+1, +1, -1};

  for (int y = 0; y < Ly_; ++y) {
    for (int x = 0; x < Lx_; ++x) {
      for (int sub = 0; sub < 2; ++sub) {
        const Sublattice s = sub == 0 ? Sublattice::A : Sublattice::B;
        const int i = site_index(x, y, s);
        for (int k = 0; k < 3; ++k) {
          int x2 = x + nnn_dx[k];
          int y2 = y + nnn_dy[k];
          int xw = 0, yw = 0, dx_wrap = 0, dy_wrap = 0;
          if (!map_coords(x2, y2, xw, yw, dx_wrap, dy_wrap)) {
            continue;
          }
          const int j = site_index(xw, yw, s);
          if (i == j) {
            continue;
          }
          const int orientation = (s == Sublattice::A) ? nnn_orient[k] : -nnn_orient[k];
          Bond b;
          b.i = i;
          b.j = j;
          b.dx_wrap = dx_wrap;
          b.dy_wrap = dy_wrap;
          b.kind = BondKind::NNN;
          b.orientation = orientation;
          nnn_bonds_.push_back(b);
        }
      }
    }
  }
}

}  // namespace ftdqmc::lattice
