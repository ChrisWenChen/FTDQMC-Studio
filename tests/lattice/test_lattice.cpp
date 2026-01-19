#include <gtest/gtest.h>
#include <unordered_map>

#include "ftdqmc/lattice/lattice.hpp"
#include "ftdqmc/lattice/chain1d.hpp"
#include "ftdqmc/lattice/square.hpp"
#include "ftdqmc/lattice/honeycomb.hpp"

static int count_dx(const std::vector<ftdqmc::Bond>& b, int dx) {
  int c = 0;
  for (auto& e : b) if (e.dx == dx) ++c;
  return c;
}

static int count_dy(const std::vector<ftdqmc::Bond>& b, int dy) {
  int c = 0;
  for (auto& e : b) if (e.dy == dy) ++c;
  return c;
}

TEST(LatticeChain1D, OBC_BondCount) {
  ftdqmc::Chain1D lat(4, false);
  EXPECT_EQ(lat.nsite(), 4);
  EXPECT_EQ(lat.nn_bonds().size(), 3u);
  EXPECT_EQ(count_dx(lat.nn_bonds(), 1), 0); // no wrap
}

TEST(LatticeChain1D, PBC_BondCountAndWrap) {
  ftdqmc::Chain1D lat(4, true);
  EXPECT_EQ(lat.nn_bonds().size(), 4u);
  EXPECT_EQ(count_dx(lat.nn_bonds(), 1), 1); // exactly one wrap bond
}

TEST(LatticeSquare, OBC_BondCount) {
  // directed bonds: +x and +y only
  // OBC => (Lx-1)*Ly + Lx*(Ly-1)
  int Lx=3, Ly=2;
  ftdqmc::SquareLattice lat(Lx, Ly, false, false);
  EXPECT_EQ(lat.nsite(), Lx*Ly);
  EXPECT_EQ(lat.nn_bonds().size(), static_cast<size_t>((Lx-1)*Ly + Lx*(Ly-1)));
  EXPECT_EQ(count_dx(lat.nn_bonds(), 1), 0);
  EXPECT_EQ(count_dy(lat.nn_bonds(), 1), 0);
}

TEST(LatticeSquare, PBC_WrapCount) {
  int Lx=3, Ly=2;
  ftdqmc::SquareLattice lat(Lx, Ly, true, true);
  // PBC directed: Lx*Ly in +x and Lx*Ly in +y => 2*Lx*Ly
  EXPECT_EQ(lat.nn_bonds().size(), static_cast<size_t>(2*Lx*Ly));
  EXPECT_EQ(count_dx(lat.nn_bonds(), 1), Ly); // for each row y, x=Lx-1 wraps once
  EXPECT_EQ(count_dy(lat.nn_bonds(), 1), Lx); // for each col x, y=Ly-1 wraps once
}

TEST(LatticeHoneycomb, SiteCountAndSublattice) {
  int Lx=2, Ly=2;
  ftdqmc::HoneycombLattice lat(Lx, Ly, true, true, false);
  EXPECT_EQ(lat.nsite(), 2*Lx*Ly);

  // construction order: (x,y,A),(x,y,B), ...
  EXPECT_EQ(lat.sublattice(0), 0);
  EXPECT_EQ(lat.sublattice(1), 1);
}

TEST(LatticeHoneycomb, NN_BondCount) {
  int Lx=2, Ly=2;
  ftdqmc::HoneycombLattice lat(Lx, Ly, true, true, false);
  // you store directed NN from A sites only, 3 bonds per cell
  EXPECT_EQ(lat.nn_bonds().size(), static_cast<size_t>(3*Lx*Ly));
}

TEST(LatticeHoneycomb, NNN_Toggle) {
  int Lx=2, Ly=2;
  ftdqmc::HoneycombLattice lat0(Lx, Ly, true, true, false);
  EXPECT_EQ(lat0.nnn_bonds().size(), 0u);

  ftdqmc::HoneycombLattice lat1(Lx, Ly, true, true, true);
  EXPECT_GT(lat1.nnn_bonds().size(), 0u);
  // with your dirs: 6 directed per site => 6 * (2*Lx*Ly)
  EXPECT_EQ(lat1.nnn_bonds().size(), static_cast<size_t>(6*2*Lx*Ly));
}
