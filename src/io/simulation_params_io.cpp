#include "ftdqmc/io/simulation_params_io.hpp"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>

namespace ftdqmc {

// Helper functions for parsing
// Trim whitespace from both ends
static inline std::string trim(std::string s) {
  auto not_space = [](unsigned char c) { return !std::isspace(c); };
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
  s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
  return s;
}

// Convert to lower case
static inline std::string to_lower(std::string s) {
  for (auto& ch : s) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
  return s;
}

// Parsing functions with error reporting
// parse_bool accepts: true/false, 1/0, yes/no, on/off (case insensitive)
static bool parse_bool(const std::string& v_raw, int line_no, const std::string& key) {
  std::string v = to_lower(trim(v_raw));
  if (v == "true" || v == "1" || v == "yes" || v == "on") return true;
  if (v == "false" || v == "0" || v == "no" || v == "off") return false;
  std::ostringstream oss;
  oss << "input.in parse error at line " << line_no << ": key '" << key
      << "' expects bool, got '" << v_raw << "'";
  throw std::runtime_error(oss.str());
}

// parse int64
static long long parse_int64(const std::string& v_raw, int line_no, const std::string& key) {
  std::string v = trim(v_raw);
  errno = 0;
  char* end = nullptr;
  long long out = std::strtoll(v.c_str(), &end, 10);
  if (errno != 0 || end == v.c_str() || *end != '\0') {
    std::ostringstream oss;
    oss << "input.in parse error at line " << line_no << ": key '" << key
        << "' expects int, got '" << v_raw << "'";
    throw std::runtime_error(oss.str());
  }
  return out;
}

// parse uint64
static unsigned long long parse_uint64(const std::string& v_raw, int line_no, const std::string& key) {
  std::string v = trim(v_raw);
  errno = 0;
  char* end = nullptr;
  unsigned long long out = std::strtoull(v.c_str(), &end, 10);
  if (errno != 0 || end == v.c_str() || *end != '\0') {
    std::ostringstream oss;
    oss << "input.in parse error at line " << line_no << ": key '" << key
        << "' expects uint64, got '" << v_raw << "'";
    throw std::runtime_error(oss.str());
  }
  return out;
}

// parse double
static double parse_double(const std::string& v_raw, int line_no, const std::string& key) {
  std::string v = trim(v_raw);
  errno = 0;
  char* end = nullptr;
  double out = std::strtod(v.c_str(), &end);
  if (errno != 0 || end == v.c_str() || *end != '\0') {
    std::ostringstream oss;
    oss << "input.in parse error at line " << line_no << ": key '" << key
        << "' expects double, got '" << v_raw << "'";
    throw std::runtime_error(oss.str());
  }
  return out;
}


// parse LatticeType
static LatticeType parse_lattice(const std::string& v_raw, int line_no) {
  std::string v = to_lower(trim(v_raw));
  if (v == "chain1d" || v == "chain" || v == "1d") return LatticeType::Chain1D;
  if (v == "square") return LatticeType::Square;
  if (v == "honeycomb" || v == "hex" || v == "hexagonal") return LatticeType::Honeycomb;

  std::ostringstream oss;
  oss << "input.in parse error at line " << line_no
      << ": invalid lattice '" << v_raw << "' (use chain1d/square/honeycomb)";
  throw std::runtime_error(oss.str());
}

// Load SimulationParams from input.in style file
SimulationParams load_simulation_params_in(const std::string& path) {
  std::ifstream fin(path);
  if (!fin) {
    throw std::runtime_error("Failed to open input file: " + path);
  }

  SimulationParams p; // defaults
  std::unordered_set<std::string> seen;

  std::string line;
  int line_no = 0;

  auto consume = [&](const std::string& key_raw, const std::string& val_raw) {
    const std::string key = to_lower(trim(key_raw));
    const std::string val = trim(val_raw);

    if (key.empty()) return;

    // Reject duplicates (helps catch mistakes)
    if (seen.count(key)) {
      std::ostringstream oss;
      oss << "input.in parse error at line " << line_no
          << ": duplicate key '" << key << "'";
      throw std::runtime_error(oss.str());
    }
    seen.insert(key);

    // Dispatch known keys
    if (key == "lattice") {
      p.lattice = parse_lattice(val, line_no);
      return;
    }
    if (key == "lx") {
      p.Lx = static_cast<int>(parse_int64(val, line_no, key));
      return;
    }
    if (key == "ly") {
      p.Ly = static_cast<int>(parse_int64(val, line_no, key));
      return;
    }
    if (key == "ltau") {
      p.Ltau = static_cast<int>(parse_int64(val, line_no, key));
      return;
    }
    if (key == "dtau") {
      p.dtau = parse_double(val, line_no, key);
      return;
    }
    if (key == "seed") {
      p.seed = static_cast<std::uint64_t>(parse_uint64(val, line_no, key));
      return;
    }
    if (key == "pbc_x") {
      p.pbc_x = parse_bool(val, line_no, key);
      return;
    }
    if (key == "pbc_y") {
      p.pbc_y = parse_bool(val, line_no, key);
      return;
    }
    if (key == "twist.enabled") {
      p.twist.enabled = parse_bool(val, line_no, key);
      return;
    }
    if (key == "twist.theta_x") {
      p.twist.theta_x = parse_double(val, line_no, key);
      return;
    }
    if (key == "twist.theta_y") {
      p.twist.theta_y = parse_double(val, line_no, key);
      return;
    }

    // Unknown key => hard error (prevents silent typos)
    std::ostringstream oss;
    oss << "input.in parse error at line " << line_no
        << ": unknown key '" << key_raw << "'";
    throw std::runtime_error(oss.str());
  };

  while (std::getline(fin, line)) {
    ++line_no;

    // Strip comments (# or ;)
    auto pos_hash = line.find('#');
    auto pos_semi = line.find(';');
    std::size_t cut = std::string::npos;
    if (pos_hash != std::string::npos) cut = pos_hash;
    if (pos_semi != std::string::npos) cut = std::min(cut, pos_semi);
    if (cut != std::string::npos) line = line.substr(0, cut);

    line = trim(line);
    if (line.empty()) continue;

    auto eq = line.find('=');
    if (eq == std::string::npos) {
      std::ostringstream oss;
      oss << "input.in parse error at line " << line_no
          << ": expected 'key = value', got '" << line << "'";
      throw std::runtime_error(oss.str());
    }

    std::string key = line.substr(0, eq);
    std::string val = line.substr(eq + 1);
    consume(key, val);
  }

  p.validate();
  return p;
}

// Helper for printing bool as "true"/"false"
static std::string bool_str(bool b) { return b ? "true" : "false"; }

// Pretty-print SimulationParams to ostream
void print_simulation_params(std::ostream& os, const SimulationParams& p) {
  // Deterministic order (good for logs + diff)
  os << "lattice = " << p.lattice_name() << "\n";
  os << "Lx = " << p.Lx << "\n";
  os << "Ly = " << p.Ly << "\n";
  os << "Ltau = " << p.Ltau << "\n";
  os << "dtau = " << p.dtau << "\n";
  os << "seed = " << static_cast<unsigned long long>(p.seed) << "\n";
  os << "pbc_x = " << bool_str(p.pbc_x) << "\n";
  os << "pbc_y = " << bool_str(p.pbc_y) << "\n";
  os << "twist.enabled = " << bool_str(p.twist.enabled) << "\n";
  os << "twist.theta_x = " << p.twist.theta_x << "\n";
  os << "twist.theta_y = " << p.twist.theta_y << "\n";
  // Useful derived values (not meant to be parsed back)
  os << "# derived: N = " << p.nsite() << ", beta = " << p.beta() << "\n";
}

// Dump SimulationParams to string
std::string dump_simulation_params_in(const SimulationParams& p) {
  std::ostringstream oss;
  print_simulation_params(oss, p);
  return oss.str();
}

// Save SimulationParams to input.in style file
void save_simulation_params_in(const SimulationParams& p, const std::string& path) {
  std::ofstream fout(path);
  if (!fout) {
    throw std::runtime_error("Failed to open output file for writing: " + path);
  }
  print_simulation_params(fout, p);
}

} // namespace ftdqmc
