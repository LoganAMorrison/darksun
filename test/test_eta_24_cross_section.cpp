/*
 * File containing tests for the 2eta -> 4eta cross section. We compare the
 * cross section with results from MadGraph
 *
 */

#include "darksun/model/boltzmann.hpp"
#include <darksun/darksun.hpp>
#include <darksun/standard_model.hpp>
#include <filesystem>
#include <fmt/core.h>
#include <fstream>
#include <gsl/gsl_errno.h>
#include <gtest/gtest.h>
#include <Eigen/Dense>

using namespace darksun;

// Conversion factor for converting pico-barn to 1/GeV^2
constexpr double PBARN_TO_GEV = 2.56819e-9;

TEST(TestModel, TestCrossSection) {
  DarkSunParameters params{10, 1e-1};
  params.mu_eta = 1.0;
  params.lec1 = 0.1;
  params.lec2 = 0.1;
  std::array<std::pair<double, double>, 5> mg_data = {
      std::make_pair(0.8 + 0.8, 0.3113 * PBARN_TO_GEV),
      std::make_pair(1.0 + 1.0, 4.331e16 * PBARN_TO_GEV),
      std::make_pair(10.0 + 10.0, 4.386e30 * PBARN_TO_GEV),
      std::make_pair(50.0 + 50.0, 2.688e40 * PBARN_TO_GEV),
      std::make_pair(100.0 + 100.0, 4.412e44 * PBARN_TO_GEV),
  };


  for (const auto &p : mg_data) {
    double cs = cross_section_2eta_4eta(p.first, params);
    fmt::print("{}, {}\n", cs, p.second);
    ASSERT_LE(std::abs(cs - p.second) / p.second, 0.1);
  }
}
