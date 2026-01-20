#include "ftdqmc/core/aux_field.hpp"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <stdexcept>

#include <hdf5.h>

namespace ftdqmc {

// --------------------------- helpers ---------------------------

static bool is_sorted_unique(const std::vector<AuxField::value_type>& v) {
  return std::is_sorted(v.begin(), v.end()) &&
         std::adjacent_find(v.begin(), v.end()) == v.end();
}

static void h5check(herr_t status, const char* msg) {
  if (status < 0) throw std::runtime_error(msg);
}

static hid_t h5open_or_create_rw(const std::string& path) {
  // Try open RW; if fails, create new.
  hid_t f = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (f >= 0) return f;
  return H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
}

static void h5write_attr_int(hid_t obj, const char* name, int value) {
  hid_t aspace = H5Screate(H5S_SCALAR);
  if (aspace < 0) throw std::runtime_error("HDF5: create scalar dataspace failed");

  // If exists, delete then recreate (simplest, robust)
  if (H5Aexists(obj, name) > 0) {
    h5check(H5Adelete(obj, name), "HDF5: delete existing attribute failed");
  }

  hid_t attr = H5Acreate2(obj, name, H5T_NATIVE_INT, aspace,
                          H5P_DEFAULT, H5P_DEFAULT);
  if (attr < 0) throw std::runtime_error("HDF5: create attribute failed");
  h5check(H5Awrite(attr, H5T_NATIVE_INT, &value), "HDF5: write attribute failed");

  H5Aclose(attr);
  H5Sclose(aspace);
}

static void h5write_attr_str(hid_t obj, const char* name, const std::string& s) {
  hid_t atype = H5Tcopy(H5T_C_S1);
  if (atype < 0) throw std::runtime_error("HDF5: copy string type failed");

  h5check(H5Tset_size(atype, s.size()), "HDF5: set string size failed");
  h5check(H5Tset_strpad(atype, H5T_STR_NULLTERM), "HDF5: set strpad failed");

  hid_t aspace = H5Screate(H5S_SCALAR);
  if (aspace < 0) throw std::runtime_error("HDF5: create scalar dataspace failed");

  if (H5Aexists(obj, name) > 0) {
    h5check(H5Adelete(obj, name), "HDF5: delete existing string attribute failed");
  }

  hid_t attr = H5Acreate2(obj, name, atype, aspace, H5P_DEFAULT, H5P_DEFAULT);
  if (attr < 0) throw std::runtime_error("HDF5: create string attribute failed");

  // Need a null-terminated buffer for HDF5 fixed-size string
  std::string tmp = s;
  tmp.push_back('\0');
  h5check(H5Awrite(attr, atype, tmp.c_str()), "HDF5: write string attribute failed");

  H5Aclose(attr);
  H5Sclose(aspace);
  H5Tclose(atype);
}

// Create group path recursively (minimal)
static hid_t h5ensure_group(hid_t file, const std::string& group_path) {
  // group_path expected like "/aux_field" or "/a/b"
  if (group_path.empty() || group_path[0] != '/') {
    throw std::invalid_argument("HDF5 group must start with '/'");
  }

  // If exists, open directly
  if (H5Lexists(file, group_path.c_str(), H5P_DEFAULT) > 0) {
    hid_t g = H5Gopen2(file, group_path.c_str(), H5P_DEFAULT);
    if (g < 0) throw std::runtime_error("HDF5: open group failed");
    return g;
  }

  // Create recursively
  hid_t current = H5Gopen2(file, "/", H5P_DEFAULT);
  if (current < 0) throw std::runtime_error("HDF5: open root group failed");

  std::string path = "/";
  std::size_t pos = 1;
  while (pos < group_path.size()) {
    std::size_t next = group_path.find('/', pos);
    std::string part = (next == std::string::npos)
                         ? group_path.substr(pos)
                         : group_path.substr(pos, next - pos);
    pos = (next == std::string::npos) ? group_path.size() : next + 1;
    if (part.empty()) continue;

    if (path.size() > 1) path += "/";
    path += part;

    if (H5Lexists(file, path.c_str(), H5P_DEFAULT) > 0) {
      hid_t g = H5Gopen2(file, path.c_str(), H5P_DEFAULT);
      H5Gclose(current);
      current = g;
    } else {
      hid_t g = H5Gcreate2(file, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(current);
      current = g;
    }
    if (current < 0) throw std::runtime_error("HDF5: create/open group component failed");
  }

  return current; // caller closes
}

// --------------------------- AuxField ---------------------------

AuxField::AuxField(int nsite, int ltau, std::vector<value_type> allowed_states)
  : nsite_(nsite), ltau_(ltau), allowed_states_(std::move(allowed_states)) {
  if (nsite_ <= 0 || ltau_ <= 0) {
    throw std::invalid_argument("AuxField: nsite and ltau must be > 0");
  }
  validate_allowed_states();
  data_.assign((std::size_t)nsite_ * (std::size_t)ltau_, (value_type)0);
}

AuxField AuxField::ising2(const SimulationParams& p) {
  return AuxField(p.nsite(), p.Ltau, std::vector<value_type>{-1, +1});
}

AuxField AuxField::ising4(const SimulationParams& p) {
  return AuxField(p.nsite(), p.Ltau, std::vector<value_type>{-2, -1, +1, +2});
}

void AuxField::validate_allowed_states() const {
  if (allowed_states_.empty()) {
    throw std::invalid_argument("AuxField: allowed_states is empty");
  }
  if (!is_sorted_unique(allowed_states_)) {
    throw std::invalid_argument("AuxField: allowed_states must be sorted & unique");
  }
  for (auto s : allowed_states_) {
    if (s == 0) {
      throw std::invalid_argument("AuxField: state 0 is not allowed (HS field typically nonzero)");
    }
  }
}

void AuxField::validate_value(value_type v) const {
  if (!std::binary_search(allowed_states_.begin(), allowed_states_.end(), v)) {
    throw std::invalid_argument("AuxField: value not in allowed_states");
  }
}

void AuxField::check_dims(int site, int tau) const {
  if (site < 0 || site >= nsite_ || tau < 0 || tau >= ltau_) {
    throw std::out_of_range("AuxField: index out of range");
  }
}

std::span<AuxField::value_type> AuxField::slice(int tau) {
  if (tau < 0 || tau >= ltau_) throw std::out_of_range("AuxField: tau out of range");
  value_type* ptr = data_.data() + (std::size_t)tau * (std::size_t)nsite_;
  return std::span<value_type>(ptr, (std::size_t)nsite_);
}

std::span<const AuxField::value_type> AuxField::slice(int tau) const {
  if (tau < 0 || tau >= ltau_) throw std::out_of_range("AuxField: tau out of range");
  const value_type* ptr = data_.data() + (std::size_t)tau * (std::size_t)nsite_;
  return std::span<const value_type>(ptr, (std::size_t)nsite_);
}

void AuxField::init(Init how,
                    const Lattice* lattice,
                    std::mt19937_64* rng,
                    value_type uniform_value,
                    value_type a_value,
                    value_type b_value) {
  switch (how) {
    case Init::Uniform:
    case Init::Ferro: {
      validate_value(uniform_value);
      std::fill(data_.begin(), data_.end(), uniform_value);
      return;
    }
    case Init::Random: {
      if (!rng) throw std::invalid_argument("AuxField::init(Random): rng is null");
      std::uniform_int_distribution<int> dist(0, (int)allowed_states_.size() - 1);
      for (auto& x : data_) x = allowed_states_[(std::size_t)dist(*rng)];
      return;
    }
    case Init::Staggered:
    case Init::AntiFerro: {
      if (!lattice) throw std::invalid_argument("AuxField::init(Staggered/AntiFerro): lattice is null");
      if (lattice->nsite() != nsite_) throw std::invalid_argument("AuxField::init: lattice.nsite mismatch");
      validate_value(a_value);
      validate_value(b_value);

      for (int tau = 0; tau < ltau_; ++tau) {
        auto row = slice(tau);
        for (int i = 0; i < nsite_; ++i) {
          row[(std::size_t)i] = (lattice->sublattice(i) == 0) ? a_value : b_value;
        }
      }
      return;
    }
  }
  throw std::runtime_error("AuxField::init: unknown Init");
}

void AuxField::init_random_from_seed(std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  init(Init::Random, nullptr, &rng);
}

// --------------------------- Binary I/O ---------------------------
// Format:
// magic "AUXFLD1\0" (8 bytes)
// u32 nsite, u32 ltau
// u32 nstates
// int8 states[nstates]
// int8 data[ltau*nsite]  (tau-major)

static void write_u32(std::ofstream& os, std::uint32_t v) {
  os.write(reinterpret_cast<const char*>(&v), sizeof(v));
}
static std::uint32_t read_u32(std::ifstream& is) {
  std::uint32_t v = 0;
  is.read(reinterpret_cast<char*>(&v), sizeof(v));
  return v;
}

void AuxField::save_binary(const std::string& path) const {
  std::ofstream os(path, std::ios::binary);
  if (!os) throw std::runtime_error("AuxField::save_binary: cannot open: " + path);

  const char magic[8] = {'A','U','X','F','L','D','1','\0'};
  os.write(magic, 8);

  write_u32(os, (std::uint32_t)nsite_);
  write_u32(os, (std::uint32_t)ltau_);
  write_u32(os, (std::uint32_t)allowed_states_.size());

  os.write(reinterpret_cast<const char*>(allowed_states_.data()),
           (std::streamsize)allowed_states_.size());

  os.write(reinterpret_cast<const char*>(data_.data()),
           (std::streamsize)data_.size());

  if (!os) throw std::runtime_error("AuxField::save_binary: write failed: " + path);
}

AuxField AuxField::load_binary(const std::string& path) {
  std::ifstream is(path, std::ios::binary);
  if (!is) throw std::runtime_error("AuxField::load_binary: cannot open: " + path);

  char magic[8];
  is.read(magic, 8);
  const char expect[8] = {'A','U','X','F','L','D','1','\0'};
  if (std::memcmp(magic, expect, 8) != 0) {
    throw std::runtime_error("AuxField::load_binary: bad magic (not AUXFLD1)");
  }

  std::uint32_t nsite = read_u32(is);
  std::uint32_t ltau  = read_u32(is);
  std::uint32_t nst   = read_u32(is);
  if (nsite == 0 || ltau == 0 || nst == 0) {
    throw std::runtime_error("AuxField::load_binary: invalid header");
  }

  std::vector<value_type> states((std::size_t)nst);
  is.read(reinterpret_cast<char*>(states.data()), (std::streamsize)states.size());

  AuxField f((int)nsite, (int)ltau, std::move(states));
  is.read(reinterpret_cast<char*>(f.data_.data()), (std::streamsize)f.data_.size());

  if (!is) throw std::runtime_error("AuxField::load_binary: read failed: " + path);
  return f;
}

// --------------------------- HDF5 I/O ---------------------------
// Layout:
// group (default "/aux_field"):
//   - dataset "s"      int8 [Ltau, Nsite]  (tau-major)
//   - dataset "states" int8 [nstates]
//   - attrs: nsite(int), ltau(int), layout("tau-major")

void AuxField::save_hdf5(const std::string& path, const std::string& group) const {
  hid_t file = h5open_or_create_rw(path);
  if (file < 0) throw std::runtime_error("AuxField::save_hdf5: cannot open/create file");

  hid_t grp = h5ensure_group(file, group);

  // If datasets exist, delete them (simplest overwrite semantics)
  if (H5Lexists(grp, "s", H5P_DEFAULT) > 0) {
    h5check(H5Ldelete(grp, "s", H5P_DEFAULT), "HDF5: delete existing dataset s failed");
  }
  if (H5Lexists(grp, "states", H5P_DEFAULT) > 0) {
    h5check(H5Ldelete(grp, "states", H5P_DEFAULT), "HDF5: delete existing dataset states failed");
  }

  // dataset "s": [ltau, nsite]
  {
    hsize_t dims[2] = {(hsize_t)ltau_, (hsize_t)nsite_};
    hid_t space = H5Screate_simple(2, dims, nullptr);
    if (space < 0) throw std::runtime_error("HDF5: create dataspace failed");

    hid_t dset = H5Dcreate2(grp, "s", H5T_NATIVE_INT8,
                            space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset < 0) throw std::runtime_error("HDF5: create dataset s failed");

    h5check(H5Dwrite(dset, H5T_NATIVE_INT8, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_.data()),
            "AuxField::save_hdf5: write dataset s failed");

    H5Dclose(dset);
    H5Sclose(space);
  }

  // dataset "states": [nstates]
  {
    hsize_t dims[1] = {(hsize_t)allowed_states_.size()};
    hid_t space = H5Screate_simple(1, dims, nullptr);
    if (space < 0) throw std::runtime_error("HDF5: create dataspace failed");

    hid_t dset = H5Dcreate2(grp, "states", H5T_NATIVE_INT8,
                            space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset < 0) throw std::runtime_error("HDF5: create dataset states failed");

    h5check(H5Dwrite(dset, H5T_NATIVE_INT8, H5S_ALL, H5S_ALL, H5P_DEFAULT, allowed_states_.data()),
            "AuxField::save_hdf5: write dataset states failed");

    H5Dclose(dset);
    H5Sclose(space);
  }

  // attributes on group
  h5write_attr_int(grp, "nsite", nsite_);
  h5write_attr_int(grp, "ltau",  ltau_);
  h5write_attr_str(grp, "layout", "tau-major");

  H5Gclose(grp);
  H5Fclose(file);
}

AuxField AuxField::load_hdf5(const std::string& path, const std::string& group) {
  hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) throw std::runtime_error("AuxField::load_hdf5: cannot open file");

  hid_t grp = H5Gopen2(file, group.c_str(), H5P_DEFAULT);
  if (grp < 0) throw std::runtime_error("AuxField::load_hdf5: cannot open group");

  // read states
  std::vector<value_type> states;
  {
    hid_t dset = H5Dopen2(grp, "states", H5P_DEFAULT);
    if (dset < 0) throw std::runtime_error("HDF5: open dataset states failed");

    hid_t space = H5Dget_space(dset);
    if (space < 0) throw std::runtime_error("HDF5: get space for states failed");

    hsize_t dims[1];
    H5Sget_simple_extent_dims(space, dims, nullptr);

    states.resize((std::size_t)dims[0]);
    h5check(H5Dread(dset, H5T_NATIVE_INT8, H5S_ALL, H5S_ALL, H5P_DEFAULT, states.data()),
            "HDF5: read states failed");

    H5Sclose(space);
    H5Dclose(dset);
  }

  // read s and infer dims
  AuxField f;
  {
    hid_t dset = H5Dopen2(grp, "s", H5P_DEFAULT);
    if (dset < 0) throw std::runtime_error("HDF5: open dataset s failed");

    hid_t space = H5Dget_space(dset);
    if (space < 0) throw std::runtime_error("HDF5: get space for s failed");

    hsize_t dims[2];
    H5Sget_simple_extent_dims(space, dims, nullptr);

    int ltau  = (int)dims[0];
    int nsite = (int)dims[1];

    f = AuxField(nsite, ltau, std::move(states));
    h5check(H5Dread(dset, H5T_NATIVE_INT8, H5S_ALL, H5S_ALL, H5P_DEFAULT, f.data_.data()),
            "HDF5: read s failed");

    H5Sclose(space);
    H5Dclose(dset);
  }

  H5Gclose(grp);
  H5Fclose(file);
  return f;
}

} // namespace ftdqmc
