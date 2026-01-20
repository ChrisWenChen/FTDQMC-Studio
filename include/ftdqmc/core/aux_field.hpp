#pragma once
#include <cstdint>
#include <random>
#include <span>
#include <string>
#include <vector>

#include "ftdqmc/core/simulation_params.hpp"
#include "ftdqmc/lattice/lattice.hpp"

namespace ftdqmc {

class AuxField {
public:
  using value_type = std::int8_t;

  // Two common auxiliary field types
  static constexpr value_type kIsing2[2] = {-1, +1};
  static constexpr value_type kIsing4[4] = {-2, -1, +1, +2};

  enum class Init : std::uint8_t {
    Uniform,
    Staggered, 
    Random,
    Ferro,
    AntiFerro
  };

  AuxField() = default;
  AuxField(int nsite, int ltau, std::vector<value_type> allowed_states);

  static AuxField ising2(const SimulationParams& p);
  static AuxField ising4(const SimulationParams& p);

  int nsite() const { return nsite_; }
  int ltau()  const { return ltau_; }
  std::size_t size() const { return data_.size(); }

  const std::vector<value_type>& allowed_states() const { return allowed_states_; }

  // tau-major: data_[tau*nsite + site]
  value_type& operator()(int site, int tau) { return data_[offset(site, tau)]; }
  value_type  operator()(int site, int tau) const { return data_[offset(site, tau)]; }

  std::span<value_type>       slice(int tau);
  std::span<const value_type> slice(int tau) const;

  std::span<value_type>       raw() { return data_; }
  std::span<const value_type> raw() const { return data_; }

  // Generation/Initialization: random, ferromagnetic, antiferromagnetic, staggered, uniform all here
  void init(Init how,
            const Lattice* lattice = nullptr,
            std::mt19937_64* rng = nullptr,
            value_type uniform_value = +1,
            value_type a_value = +1,
            value_type b_value = -1);

  void init_random_from_seed(std::uint64_t seed);

  // Save/Load (initially zero-dependency binary, later can switch to HDF5 without changing interface)
  void save_binary(const std::string& path) const;
  static AuxField load_binary(const std::string& path);

private:
  int nsite_ = 0;
  int ltau_  = 0;
  std::vector<value_type> allowed_states_;
  std::vector<value_type> data_; // [tau][site]

  std::size_t offset(int site, int tau) const {
    return (std::size_t)tau * (std::size_t)nsite_ + (std::size_t)site;
  }

  void validate_allowed_states() const;
  void validate_value(value_type v) const;
};

} // namespace ftdqmc
