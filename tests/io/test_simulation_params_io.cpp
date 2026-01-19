#include <gtest/gtest.h>

#include <fstream>
#include <string>

#include "ftdqmc/io/simulation_params_io.hpp"

static std::string write_temp_in(const std::string& content) {
  const std::string path = "test_input_params.in";
  std::ofstream out(path, std::ios::binary);
  out << content;
  out.close();
  return path;
}

TEST(ParamsIO, LoadAndDump) {
  const char* text = R"(
# Example input.in
lattice = chain1d
Lx = 12
Ly = 1

Ltau = 200
dtau = 0.05
seed = 123456

pbc_x = true
pbc_y = false

twist.enabled = true
twist.theta_x = 0.3
twist.theta_y = 0.0
)";

  const std::string path = write_temp_in(text);
  auto p = ftdqmc::load_simulation_params_in(path);

  EXPECT_EQ(p.lattice_name(), "chain1d");
  EXPECT_EQ(p.Lx, 12);
  EXPECT_EQ(p.Ly, 1);
  EXPECT_EQ(p.Ltau, 200);
  EXPECT_DOUBLE_EQ(p.dtau, 0.05);
  EXPECT_EQ(p.seed, 123456ULL);
  EXPECT_TRUE(p.pbc_x);
  EXPECT_FALSE(p.pbc_y);
  EXPECT_TRUE(p.twist.enabled);
  EXPECT_DOUBLE_EQ(p.twist.theta_x, 0.3);
  EXPECT_DOUBLE_EQ(p.twist.theta_y, 0.0);

  const std::string dumped = ftdqmc::dump_simulation_params_in(p);
  EXPECT_NE(dumped.find("lattice = chain1d"), std::string::npos);
  EXPECT_NE(dumped.find("Lx = 12"), std::string::npos);
  EXPECT_NE(dumped.find("twist.enabled = true"), std::string::npos);
  EXPECT_NE(dumped.find("# derived:"), std::string::npos);
}

TEST(ParamsIO, SaveAndReloadRoundtrip) {
  ftdqmc::SimulationParams p;
  p.lattice = ftdqmc::LatticeType::Honeycomb;
  p.Lx = 4;
  p.Ly = 3;
  p.Ltau = 40;
  p.dtau = 0.1;
  p.seed = 999;
  p.pbc_x = true;
  p.pbc_y = true;
  p.twist.enabled = false;
  p.twist.theta_x = 0.0;
  p.twist.theta_y = 0.0;
  p.validate();

  const std::string outpath = "test_saved_params.in";
  ftdqmc::save_simulation_params_in(p, outpath);

  auto q = ftdqmc::load_simulation_params_in(outpath);

  EXPECT_EQ(q.lattice_name(), p.lattice_name());
  EXPECT_EQ(q.Lx, p.Lx);
  EXPECT_EQ(q.Ly, p.Ly);
  EXPECT_EQ(q.Ltau, p.Ltau);
  EXPECT_DOUBLE_EQ(q.dtau, p.dtau);
  EXPECT_EQ(q.seed, p.seed);
  EXPECT_EQ(q.pbc_x, p.pbc_x);
  EXPECT_EQ(q.pbc_y, p.pbc_y);
  EXPECT_EQ(q.twist.enabled, p.twist.enabled);
  EXPECT_DOUBLE_EQ(q.twist.theta_x, p.twist.theta_x);
  EXPECT_DOUBLE_EQ(q.twist.theta_y, p.twist.theta_y);
}
