// Simple script to output data on the 2eta -> 4eta thermal cross section

#include <cmath>
#include <darksun/darksun.hpp>
#include <filesystem>
#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fstream>

using namespace darksun;

int main() {
  DarkSunParameters params{10, 1e-1};

  constexpr size_t num_xs = 200;
  constexpr double logx_min = -1.0;
  constexpr double logx_max = 2.0;
  constexpr double logx_stp = (logx_max - logx_min) / double(num_xs - 1);

  const std::string file =
      std::filesystem::current_path().append("../rundata/tc_data.csv");
  std::ofstream ofile(file);

  ofile << "X,TC24,TC42\n";

  const double coeff = pow(256.0 * pow(M_PI, 4) * pow(m_eta(params), 7) /
                               (pow(params.lam, 8) * pow(params.n, 2)),
                           2);

  for (size_t i = 0; i < num_xs; i++) {
    const double logx = logx_min + i * logx_stp;
    const double x = pow(10.0, logx);
    const double tc24 = thermal_cross_section_2eta_4eta(x, params) / coeff;
    const double tc42 = thermal_cross_section_4eta_2eta(x, params) / coeff;

    fmt::print(ofile, "{},{},{}\n", x, tc24, tc42);
  }
  ofile.close();
}
