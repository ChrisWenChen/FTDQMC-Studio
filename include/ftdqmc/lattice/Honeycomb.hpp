#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace ftdqmc::lattice {

// Site indexing (fixed):
// site = 2 * (x + Lx * y) + sub, sub=0(A), 1(B)
// Ns = 2 * Lx * Ly

enum class BoundaryCondition { OBC, PBC, TBC };

enum class Sublattice { A = 0, B = 1 };

enum class BondKind { NN, NNN };

struct Bond {
  int i = 0;
  int j = 0;
  int dx_wrap = 0;  // wrap count in x (for PBC/TBC)
  int dy_wrap = 0;  // wrap count in y (for PBC/TBC)
  BondKind kind = BondKind::NN;
  int orientation = 0;  // For Haldane NNN: +1 or -1; 0 for NN
};

class HoneycombLattice {
 public:
  HoneycombLattice(int Lx, int Ly, BoundaryCondition bc,
                   double twist_theta_x = 0.0, double twist_theta_y = 0.0);

  int Ns() const { return ns_; }
  int Lx() const { return Lx_; }
  int Ly() const { return Ly_; }

  BoundaryCondition bc() const { return bc_; }
  double twist_theta_x() const { return twist_theta_x_; }
  double twist_theta_y() const { return twist_theta_y_; }

  Sublattice sublattice(int site) const;

  const std::vector<Bond>& nn_bonds() const { return nn_bonds_; }
  const std::vector<Bond>& nnn_bonds() const { return nnn_bonds_; }

 private:
  int Lx_ = 0;
  int Ly_ = 0;
  int ns_ = 0;
  BoundaryCondition bc_ = BoundaryCondition::OBC;
  double twist_theta_x_ = 0.0;
  double twist_theta_y_ = 0.0;

  std::vector<Bond> nn_bonds_;
  std::vector<Bond> nnn_bonds_;

  int site_index(int x, int y, Sublattice sub) const;
  bool map_coords(int x, int y, int& xw, int& yw, int& dx_wrap, int& dy_wrap) const;
  void build_bonds();
};

}  // namespace ftdqmc::lattice
