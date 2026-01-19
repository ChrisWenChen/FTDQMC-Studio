#pragma once
#include <vector>

namespace ftdqmc {

// Bond: i -> j, and (dx,dy) indicates wrapping across boundary in unit-cell steps.
// This is crucial for twisted boundary condition later.
struct Bond {
  int i;
  int j;
  int dx;
  int dy;
};

struct Site {
  int x;
  int y;
  int sub; // sublattice index: chain/square=0, honeycomb A=0/B=1
};

class Lattice {
public:
  virtual ~Lattice() = default;

  int nsite() const { return static_cast<int>(sites_.size()); }
  const std::vector<Site>& sites() const { return sites_; }

  const std::vector<Bond>& nn_bonds() const { return nn_; }   // nearest neighbors
  const std::vector<Bond>& nnn_bonds() const { return nnn_; } // next-nearest neighbors (Haldane uses this)

  int sublattice(int i) const { return sites_[i].sub; }

protected:
  std::vector<Site> sites_;
  std::vector<Bond> nn_;
  std::vector<Bond> nnn_;
};

} // namespace ftdqmc
