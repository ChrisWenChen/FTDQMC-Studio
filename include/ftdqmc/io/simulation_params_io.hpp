#pragma once
#include <iosfwd>
#include <string>

#include "ftdqmc/core/simulation_params.hpp"

namespace ftdqmc {

// Load from a key=value input file (input.in style)
SimulationParams load_simulation_params_in(const std::string& path);

// Save to a key=value input file (deterministic order)
void save_simulation_params_in(const SimulationParams& p, const std::string& path);

// Dump to text (same format as save, but returns a string)
std::string dump_simulation_params_in(const SimulationParams& p);

// Convenience printing
void print_simulation_params(std::ostream& os, const SimulationParams& p);

} // namespace ftdqmc
