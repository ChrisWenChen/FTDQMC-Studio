#pragma once
#include <cstdint>
#include <string>
#include <stdexcept>

namespace ftdqmc {

enum class LatticeType { Chain1D, Square, Honeycomb };


// Twisted boundary conditions (TBC):
// theta_x/y are in radians. Default 0 => ordinary PBC.
// In practice, phase is applied to hoppings that wrap around boundaries.
struct TwistBC {
  bool enabled = false;
  double theta_x = 0.0;
  double theta_y = 0.0;
};

struct SimulationParams {
  // Lattice
  LatticeType lattice = LatticeType::Honeycomb;
  int Lx = 4;
  int Ly = 4;

  // Imaginary time discretization
  int Ltau = 40;      // number of time slices
  double dtau = 0.1;  // trotter step
  double beta() const { return Ltau * dtau; }

  // RNG
  std::uint64_t seed = 12345;

  // Boundary conditions (future: twist angles)
  bool pbc_x = true;
  bool pbc_y = true;

  // Twisted boundary condition (optional)
  TwistBC twist;


  // Derived: number of unit cells
  int ncell() const {
    switch (lattice) {
      case LatticeType::Chain1D:  return Lx;
      case LatticeType::Square:   return Lx * Ly;
      case LatticeType::Honeycomb:return Lx * Ly;
    }
    return 0;
  }

  // Derived: number of sites (spinless, single-orbital)
  // chain: 1 site / cell
  // square: 1 site / cell
  // honeycomb: 2 sites / cell (A/B)
  int nsite() const {
    switch (lattice) {
      case LatticeType::Chain1D:  return ncell();
      case LatticeType::Square:   return ncell();
      case LatticeType::Honeycomb:return 2 * ncell();
    }
    return 0;
  }

  void validate() const {
    if (Lx <= 0) throw std::invalid_argument("Lx must be positive.");
    if (lattice != LatticeType::Chain1D && Ly <= 0)
      throw std::invalid_argument("Ly must be positive for 2D lattices.");
    if (lattice == LatticeType::Chain1D && Ly != 1)
      throw std::invalid_argument("For Chain1D, set Ly=1 (recommended).");

    if (Ltau <= 0) throw std::invalid_argument("Ltau must be positive.");
    if (dtau <= 0.0) throw std::invalid_argument("dtau must be positive.");

    // Twist sanity: if twist enabled, usually you want PBC in that direction
    if (twist.enabled) {
      if (!pbc_x && twist.theta_x != 0.0)
        throw std::invalid_argument("TwistBC: theta_x != 0 requires pbc_x=true.");
      if (!pbc_y && twist.theta_y != 0.0)
        throw std::invalid_argument("TwistBC: theta_y != 0 requires pbc_y=true.");
    }
  }

  std::string lattice_name() const {
    switch (lattice) {
      case LatticeType::Chain1D:  return "chain1d";
      case LatticeType::Square:   return "square";
      case LatticeType::Honeycomb:return "honeycomb";
    }
    return "unknown";
  }
};

} // namespace ftdqmc