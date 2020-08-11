//
// Created by logan on 8/9/20.
//

#include <darksun/darksun.hpp>
#include <darksun/standard_model.hpp>
#include <filesystem>
#include <fmt/core.h>
#include <fstream>
#include <gsl/gsl_errno.h>
#include <gtest/gtest.h>

using namespace darksun;

TEST(TestModel, TestCrossSection) {
  DarkSunParameters params{10, 1e-1};
  params.mu_eta = 1.0;
  params.lec1 = 0.1;
  params.lec2 = 0.1;
  std::array<std::pair<double, double>, 5> mg_data = {
      std::make_pair(1.0 + 1.0, 1.10971e8),
      std::make_pair(5.0 + 5.0, 6.88532e17),
      std::make_pair(10.0 + 10.0, 1.12898e22),
      std::make_pair(50.0 + 50.0, 6.89559e31),
      std::make_pair(100.0 + 100.0, 1.1289836e36),
  };

  for (const auto &p : mg_data) {
    double cs = cross_section_2eta_4eta(p.first, params);
    fmt::print("{}, {}\n", cs, p.second);
    ASSERT_LE(std::abs(cs - p.second) / p.second, 0.1);
  }
}

TEST(TestModel, TestThermalCrossSection) {
  gsl_set_error_handler_off();
  std::string file =
      std::filesystem::current_path().append("../rundata/tc_data.csv");
  std::ofstream ofile(file);
  DarkSunParameters params{10, 1e-1};
  const size_t num_xs = 150.0;
  const double logxmin = -1.0;
  const double logxmax = 2.0;
  const double logxstp = (logxmax - logxmin) / double(logxstp - 1);

  ofile << "X,TC24,TC42\n";

  for (size_t i = 0; i < num_xs; i++) {
    const double logx = logxmin + i * logxstp;
    const double x = pow(10.0, logx);
    const double tc24 = thermal_cross_section_2eta_4eta(x, params);
    const double tc42 = thermal_cross_section_4eta_2eta(x, params);
    ofile << x << "," << tc24 << "," << tc42 << "\n";
  }
  ofile.close();
}

TEST(TestModel, TestEquilibrium) {
  DarkSunParameters params{10, 1e-1};

  const double td = 1.0;

  const double heff = dark_heff(td, params);
  const double geff = dark_geff(td, params);
  const double neq_e = neq_eta(td, params);
  const double neq_d = neq_del(td, params);
  const double yeq_e = yeq_eta(td, 1.0, params);
  const double yeq_d = yeq_del(td, 1.0, params);
  const double weq_e = weq_eta(td, 1.0, params);
  const double weq_d = weq_del(td, 1.0, params);

  fmt::print("td = {}, heff = {}, geff = {}\n\n", td, heff, geff);

  fmt::print("td = {}, neqe = {}, neqe/s = {}, log(neqe/s) = {}\n", td, neq_e,
             neq_e / StandardModel::entropy_density(td),
             log(neq_e / StandardModel::entropy_density(td)));
  fmt::print("yeqe = {}, weq_e = {}\n\n", yeq_e, weq_e);

  fmt::print("td = {}, neqd = {}, neqd/s = {}, log(neqd/s) = {}\n", td, neq_e,
             neq_d / StandardModel::entropy_density(td),
             log(neq_d / StandardModel::entropy_density(td)));
  fmt::print("yeqd = {}, weq_d = {}\n\n", yeq_d, weq_d);
}

TEST(TestModel, TestSigma24) {
  DarkSunParameters params{10, 1e-1};
  fmt::print("sigma_24(cme=1GeV) = {}\n",
             cross_section_2eta_4eta(20.0, params));

  fmt::print("sigma_24(cme=1GeV) = {}\n",
             cross_section_2eta_4eta(99.0 * m_eta(params), params));
  fmt::print("sigma_24(cme=1GeV) = {}\n",
             cross_section_2eta_4eta(101.0 * m_eta(params), params));
}

TEST(TestModel, TestSigmav) {
  DarkSunParameters params{10, 1e-1};
  fmt::print("<sigma*v>_24(x=10) = {}\n",
             thermal_cross_section_2eta_4eta(10, params));
  fmt::print("<sigma*v>_42(x=10) = {}\n",
             thermal_cross_section_4eta_2eta(10, params));

  fmt::print("<sigma*v>_ed(x=1) = {}\n",
             thermal_cross_section_2eta_2del(1, params));
  fmt::print("<sigma*v>_ed(x=5) = {}\n",
             thermal_cross_section_2eta_2del(5, params));
  fmt::print("<sigma*v>_ed(x=10) = {}\n",
             thermal_cross_section_2eta_2del(10, params));
  fmt::print("<sigma*v>_ed(x=50) = {}\n",
             thermal_cross_section_2eta_2del(50, params));
}

TEST(TestModel, TestSolveBoltzmann) {
  DarkSunParameters params{10, 1e-1};

  solve_boltzmann(1e-9, 1e-9, params);

  for (size_t i = 0; i < params.ts.size(); i++) {
    fmt::print("logx = {}, weta = {}, ydel = {}\n", params.ts[i],
               params.ys[i][0], params.ys[i][1]);
  }
}
