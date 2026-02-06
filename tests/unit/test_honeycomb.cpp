#include <gtest/gtest.h>

#include <algorithm>
#include <set>
#include <vector>

#include "ftdqmc/lattice/Honeycomb.hpp"
#include "ftdqmc/linalg/DenseMatrix.hpp"
#include "ftdqmc/model/KBuilder.hpp"
#include "ftdqmc/model/Model.hpp"

namespace {
void check_sublattice_consistency(const ftdqmc::lattice::HoneycombLattice& lat) {
  int ns = lat.Ns();
  for (int site = 0; site < ns; ++site) {
    auto sub = lat.sublattice(site);
    if (site % 2 == 0) {
      EXPECT_EQ(sub, ftdqmc::lattice::Sublattice::A);
    } else {
      EXPECT_EQ(sub, ftdqmc::lattice::Sublattice::B);
    }
  }
}

void check_coordination_numbers(const ftdqmc::lattice::HoneycombLattice& lat) {
  int ns = lat.Ns();
  std::vector<int> nn_count(ns, 0);
  std::vector<int> nnn_count(ns, 0);

  for (const auto& bond : lat.nn_bonds()) {
    nn_count[bond.i]++;
    nn_count[bond.j]++;
  }

  for (const auto& bond : lat.nnn_bonds()) {
    nnn_count[bond.i]++;
    nnn_count[bond.j]++;
  }

  for (int site = 0; site < ns; ++site) {
    if (lat.bc() == ftdqmc::lattice::BoundaryCondition::PBC) {
      EXPECT_EQ(nn_count[site], 3);
      EXPECT_EQ(nnn_count[site], 6);
    }
  }
}

void check_no_duplicate_bonds(const ftdqmc::lattice::HoneycombLattice& lat) {
  std::set<std::pair<int, int>> nn_pairs;
  std::set<std::pair<int, int>> nnn_pairs;

  for (const auto& bond : lat.nn_bonds()) {
    int i = std::min(bond.i, bond.j);
    int j = std::max(bond.i, bond.j);
    auto pair = std::make_pair(i, j);
    EXPECT_TRUE(nn_pairs.find(pair) == nn_pairs.end());
    nn_pairs.insert(pair);
  }

  for (const auto& bond : lat.nnn_bonds()) {
    int i = std::min(bond.i, bond.j);
    int j = std::max(bond.i, bond.j);
    auto pair = std::make_pair(i, j);
    EXPECT_TRUE(nnn_pairs.find(pair) == nnn_pairs.end());
    nnn_pairs.insert(pair);
  }
}

std::vector<int> collect_nn(const ftdqmc::lattice::HoneycombLattice& lat, int site) {
  std::vector<int> neighbors;
  for (const auto& bond : lat.nn_bonds()) {
    if (bond.i == site) neighbors.push_back(bond.j);
    if (bond.j == site) neighbors.push_back(bond.i);
  }
  std::sort(neighbors.begin(), neighbors.end());
  return neighbors;
}

std::vector<int> collect_nnn(const ftdqmc::lattice::HoneycombLattice& lat, int site) {
  std::vector<int> neighbors;
  for (const auto& bond : lat.nnn_bonds()) {
    if (bond.i == site) neighbors.push_back(bond.j);
    if (bond.j == site) neighbors.push_back(bond.i);
  }
  std::sort(neighbors.begin(), neighbors.end());
  return neighbors;
}

void check_3x3_neighbors(const ftdqmc::lattice::HoneycombLattice& lat) {
  EXPECT_EQ(lat.Ns(), 18);
  for (int site = 0; site < lat.Ns(); ++site) {
    auto nn = collect_nn(lat, site);
    auto nnn = collect_nnn(lat, site);
    EXPECT_EQ(nn.size(), 3u);
    EXPECT_EQ(nnn.size(), 6u);
  }
}

struct Coord {
  int x = 0;
  int y = 0;
  ftdqmc::lattice::Sublattice sub = ftdqmc::lattice::Sublattice::A;
};

Coord site_to_coord(int site, int Lx) {
  Coord c;
  c.sub = (site % 2 == 0) ? ftdqmc::lattice::Sublattice::A : ftdqmc::lattice::Sublattice::B;
  int cell = site / 2;
  c.x = cell % Lx;
  c.y = cell / Lx;
  return c;
}

int expected_nnn_orientation(const Coord& from, const Coord& to, int Lx, int Ly, int dx_wrap, int dy_wrap) {
  const int dx = (to.x - from.x) + dx_wrap * Lx;
  const int dy = (to.y - from.y) + dy_wrap * Ly;

  // Allowed NNN displacements in builder: (1,0), (0,-1), (-1,1)
  int base = 0;
  if (dx == 1 && dy == 0) {
    base = +1;
  } else if (dx == 0 && dy == -1) {
    base = +1;
  } else if (dx == -1 && dy == 1) {
    base = -1;
  } else {
    return 0;  // not a valid NNN displacement for this orientation convention
  }

  return (from.sub == ftdqmc::lattice::Sublattice::A) ? base : -base;
}

int site_index_from_coord(int x, int y, int Lx, int Ly, ftdqmc::lattice::Sublattice sub) {
  int xw = x % Lx;
  int yw = y % Ly;
  if (xw < 0) xw += Lx;
  if (yw < 0) yw += Ly;
  return 2 * (xw + Lx * yw) + static_cast<int>(sub);
}
}

TEST(HoneycombLattice, NsMatchesFormula) {
  ftdqmc::lattice::HoneycombLattice lat(2, 3, ftdqmc::lattice::BoundaryCondition::PBC);
  EXPECT_EQ(lat.Ns(), 2 * lat.Lx() * lat.Ly());
  EXPECT_EQ(lat.Ns(), 12);
}

TEST(HoneycombLattice, NNBondCountSmallPBC) {
  ftdqmc::lattice::HoneycombLattice lat(1, 1, ftdqmc::lattice::BoundaryCondition::PBC);
  EXPECT_EQ(static_cast<int>(lat.nn_bonds().size()), 3);
}

TEST(HoneycombLattice, NNBondCountTwoByTwoPBC) {
  ftdqmc::lattice::HoneycombLattice lat(2, 2, ftdqmc::lattice::BoundaryCondition::PBC);
  int expected = 3 * 2 * 2;
  EXPECT_EQ(static_cast<int>(lat.nn_bonds().size()), expected);
}

TEST(HoneycombLattice, NNNBondCountTwoByTwoPBC) {
  ftdqmc::lattice::HoneycombLattice lat(2, 2, ftdqmc::lattice::BoundaryCondition::PBC);
  int expected = 3 * lat.Ns();
  EXPECT_EQ(static_cast<int>(lat.nnn_bonds().size()), expected);
}

TEST(HoneycombLattice, SublatticeConsistency) {
  ftdqmc::lattice::HoneycombLattice lat(3, 3, ftdqmc::lattice::BoundaryCondition::PBC);
  check_sublattice_consistency(lat);
}

TEST(HoneycombLattice, CoordinationNumbersPBC) {
  ftdqmc::lattice::HoneycombLattice lat(3, 3, ftdqmc::lattice::BoundaryCondition::PBC);
  check_coordination_numbers(lat);
}

TEST(HoneycombLattice, NoDuplicateBondsPBC) {
  ftdqmc::lattice::HoneycombLattice lat(4, 4, ftdqmc::lattice::BoundaryCondition::PBC);
  check_no_duplicate_bonds(lat);
}

TEST(HoneycombLattice, OBCHasFewerBonds) {
  ftdqmc::lattice::HoneycombLattice lat_obc(2, 2, ftdqmc::lattice::BoundaryCondition::OBC);
  ftdqmc::lattice::HoneycombLattice lat_pbc(2, 2, ftdqmc::lattice::BoundaryCondition::PBC);
  EXPECT_LT(lat_obc.nn_bonds().size(), lat_pbc.nn_bonds().size());
}

TEST(HoneycombLattice, MatrixDimensionsMatchSites) {
  using namespace ftdqmc;
  lattice::HoneycombLattice lat(3, 3, lattice::BoundaryCondition::PBC);
  model::ModelSpec spec;
  spec.type = model::ModelType::Hubbard;
  spec.hubbard.t1 = 1.0;
  spec.hubbard.mu = 0.0;
  linalg::DenseMatrix K = model::build_dense(lat, spec);
  EXPECT_EQ(K.rows(), lat.Ns());
  EXPECT_EQ(K.cols(), lat.Ns());
}

TEST(HoneycombLattice, DifferentLatticeSizes) {
  for (int Lx = 1; Lx <= 4; ++Lx) {
    for (int Ly = 1; Ly <= 4; ++Ly) {
      ftdqmc::lattice::HoneycombLattice lat(Lx, Ly, ftdqmc::lattice::BoundaryCondition::PBC);
      EXPECT_EQ(lat.Ns(), 2 * Lx * Ly);
      EXPECT_EQ(lat.Lx(), Lx);
      EXPECT_EQ(lat.Ly(), Ly);
    }
  }
}

TEST(HoneycombLattice, BondWrappingInPBC) {
  ftdqmc::lattice::HoneycombLattice lat(2, 2, ftdqmc::lattice::BoundaryCondition::PBC);
  bool has_wrapped = false;
  for (const auto& bond : lat.nn_bonds()) {
    if (bond.dx_wrap != 0 || bond.dy_wrap != 0) {
      has_wrapped = true;
      break;
    }
  }
  EXPECT_TRUE(has_wrapped);
}

TEST(HoneycombLattice, BondKindCorrectness) {
  ftdqmc::lattice::HoneycombLattice lat(2, 2, ftdqmc::lattice::BoundaryCondition::PBC);
  for (const auto& bond : lat.nn_bonds()) {
    EXPECT_EQ(bond.kind, ftdqmc::lattice::BondKind::NN);
  }
  for (const auto& bond : lat.nnn_bonds()) {
    EXPECT_EQ(bond.kind, ftdqmc::lattice::BondKind::NNN);
  }
}

TEST(HoneycombLattice, Hardcoded3x3NeighborsPBC) {
  ftdqmc::lattice::HoneycombLattice lat(3, 3, ftdqmc::lattice::BoundaryCondition::PBC);

  struct Case { int site; std::vector<int> expected; };
  std::vector<Case> cases = {
    {0, {1, 5, 13}},
    {1, {0, 2, 6}},
    {2, {1, 3, 15}},
    {8, {3, 7, 9}},
    {17, {4, 12, 16}},
  };

  for (const auto& tc : cases) {
    auto expected = tc.expected;
    std::sort(expected.begin(), expected.end());
    auto actual = collect_nn(lat, tc.site);
    ASSERT_EQ(actual.size(), expected.size());
    EXPECT_EQ(actual, expected);
  }

  check_3x3_neighbors(lat);
}

TEST(HoneycombLattice, NNNPhaseOrientationPBC3x3) {
  ftdqmc::lattice::HoneycombLattice lat(3, 3, ftdqmc::lattice::BoundaryCondition::PBC);

  for (const auto& bond : lat.nnn_bonds()) {
    const auto ci = site_to_coord(bond.i, lat.Lx());
    const auto cj = site_to_coord(bond.j, lat.Lx());

    const int expected = expected_nnn_orientation(ci, cj, lat.Lx(), lat.Ly(),
                                                   bond.dx_wrap, bond.dy_wrap);
    ASSERT_NE(expected, 0);
    EXPECT_EQ(bond.orientation, expected);
  }
}

TEST(HoneycombLattice, Hardcoded3x3NNNForSite0) {
  ftdqmc::lattice::HoneycombLattice lat(3, 3, ftdqmc::lattice::BoundaryCondition::PBC);

  const int site = 0;
  const auto nnn = collect_nnn(lat, site);
  ASSERT_EQ(nnn.size(), 6u);

  // Hardcoded NNN neighbors for site 0 (A sublattice) in 3x3 PBC.
  std::vector<int> expected = {2, 4, 6, 10, 12, 14};
  std::sort(expected.begin(), expected.end());
  auto actual = nnn;
  std::sort(actual.begin(), actual.end());
  EXPECT_EQ(actual, expected);

  // Check orientation (phase sign) for each NNN bond from site 0.
  for (const auto& bond : lat.nnn_bonds()) {
    if (bond.i != site) {
      continue;
    }
    int other = bond.j;
    // Determine expected orientation using displacement rules.
    const auto ci = site_to_coord(site, lat.Lx());
    const auto cj = site_to_coord(other, lat.Lx());
    const int expected_orient = expected_nnn_orientation(ci, cj, lat.Lx(), lat.Ly(),
                                                         bond.dx_wrap, bond.dy_wrap);
    ASSERT_NE(expected_orient, 0);
    EXPECT_EQ(bond.orientation, expected_orient);
  }
}

TEST(HoneycombLattice, NeighborSymmetry) {
  ftdqmc::lattice::HoneycombLattice lat(3, 3, ftdqmc::lattice::BoundaryCondition::PBC);

  for (int site = 0; site < lat.Ns(); ++site) {
    auto nn = collect_nn(lat, site);
    for (int neighbor : nn) {
      auto nn_of_neighbor = collect_nn(lat, neighbor);
      EXPECT_NE(std::find(nn_of_neighbor.begin(), nn_of_neighbor.end(), site), nn_of_neighbor.end());
    }
  }

  for (int site = 0; site < lat.Ns(); ++site) {
    auto nnn = collect_nnn(lat, site);
    for (int neighbor : nnn) {
      auto nnn_of_neighbor = collect_nnn(lat, neighbor);
      EXPECT_NE(std::find(nnn_of_neighbor.begin(), nnn_of_neighbor.end(), site), nnn_of_neighbor.end());
    }
  }
}

TEST(HoneycombLattice, NNBondsCrossSublattices) {
  ftdqmc::lattice::HoneycombLattice lat(3, 3, ftdqmc::lattice::BoundaryCondition::PBC);
  for (const auto& bond : lat.nn_bonds()) {
    auto sub_i = lat.sublattice(bond.i);
    auto sub_j = lat.sublattice(bond.j);
    EXPECT_NE(sub_i, sub_j);
  }
}

TEST(HoneycombLattice, NNNBondsSameSublattice) {
  ftdqmc::lattice::HoneycombLattice lat(3, 3, ftdqmc::lattice::BoundaryCondition::PBC);
  for (const auto& bond : lat.nnn_bonds()) {
    auto sub_i = lat.sublattice(bond.i);
    auto sub_j = lat.sublattice(bond.j);
    EXPECT_EQ(sub_i, sub_j);
  }
}

TEST(KBuilder, HermitianHubbardAndHaldane) {
  using namespace ftdqmc;
  lattice::HoneycombLattice lat(1, 1, lattice::BoundaryCondition::PBC);

  model::ModelSpec hubbard;
  hubbard.type = model::ModelType::Hubbard;
  hubbard.hubbard.t1 = 1.0;
  hubbard.hubbard.mu = 0.2;
  linalg::DenseMatrix K1 = model::build_dense(lat, hubbard);
  double err1 = linalg::hermitian_error_norm(K1);
  EXPECT_LT(err1, 1e-12);

  model::ModelSpec haldane;
  haldane.type = model::ModelType::HaldaneHubbard;
  haldane.haldane.t1 = 1.0;
  haldane.haldane.t2 = 0.15;
  haldane.haldane.phi = 1.234;
  haldane.haldane.mu = -0.1;
  haldane.haldane.staggered_mass = 0.05;
  linalg::DenseMatrix K2 = model::build_dense(lat, haldane);
  double err2 = linalg::hermitian_error_norm(K2);
  EXPECT_LT(err2, 1e-12);
}
