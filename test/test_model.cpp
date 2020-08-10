//
// Created by logan on 8/9/20.
//

#include <darksun/darksun.hpp>
#include <darksun/standard_model.hpp>
#include <fmt/core.h>
#include <gtest/gtest.h>

using namespace darksun::model;
using namespace darksun;

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
