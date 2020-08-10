
#include <boost/timer/progress_display.hpp>
#include <darksun/phase_space.hpp>
#include <filesystem>
#include <fmt/core.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>

using namespace darksun;

double amp_4pt(const FourMomentum &, const FourMomentum &, const FourMomentum &,
               const FourMomentum &, const FourMomentum &,
               const FourMomentum &);

double amp_6pt(const FourMomentum &, const FourMomentum &, const FourMomentum &,
               const FourMomentum &, const FourMomentum &,
               const FourMomentum &);

std::tuple<double, double, double> scaled_cross_section(double, size_t);

int main() {
  std::string fname_zs =
      std::filesystem::current_path().append("../rundata/cs_data/log10_zs.dat");
  std::string fname_cs44 = std::filesystem::current_path().append(
      "../rundata/cs_data/log10_cs44.dat");
  std::string fname_cs66 = std::filesystem::current_path().append(
      "../rundata/cs_data/log10_cs66.dat");
  std::string fname_cs46 = std::filesystem::current_path().append(
      "../rundata/cs_data/log10_cs46.dat");

  std::ofstream file_zs;
  std::ofstream file_cs44;
  std::ofstream file_cs66;
  std::ofstream file_cs46;
  file_zs.open(fname_zs);
  file_cs44.open(fname_cs44);
  file_cs66.open(fname_cs66);
  file_cs46.open(fname_cs46);

  file_zs << std::setprecision(17);
  file_cs44 << std::setprecision(17);
  file_cs66 << std::setprecision(17);
  file_cs46 << std::setprecision(17);

  const double logz_min = log10(4.0 + 1e-5);
  const double logz_max = log10(100.0);
  const size_t num_zs = 500;
  const double logz_step = (logz_max - logz_min) / double(num_zs - 1);
  boost::timer::progress_display progress(num_zs);

  for (size_t i = 0; i < num_zs; i++) {
    const double logz = logz_min + i * logz_step;
    const double z = pow(10.0, logz);
    const auto res = scaled_cross_section(z, 1000000);
    file_zs << logz << "\n";
    file_cs44 << log10(std::get<0>(res)) << "\n";
    file_cs66 << log10(std::get<1>(res)) << "\n";
    file_cs46 << log10(std::get<2>(res)) << "\n";
    ++progress;
  }

  file_zs.close();
  file_cs44.close();
  file_cs66.close();
  file_cs46.close();

  return 0;
}

//===========================================================================
//---- Amplitude for 2eta->4eta using only 4pt interactions -----------------
//===========================================================================

double amp_4pt(const FourMomentum &p1, const FourMomentum &p2,
               const FourMomentum &p3, const FourMomentum &p4,
               const FourMomentum &p5, const FourMomentum &p6) {
  return ((scalar_product(p2, p2) * scalar_product(p3, p5) -
           scalar_product(p2, p5) *
               (scalar_product(p3, p3) + 2 * scalar_product(p3, p5)) +
           scalar_product(p2, p3) *
               (2 * scalar_product(p2, p5) - 2 * scalar_product(p3, p5) -
                scalar_product(p5, p5))) *
          (scalar_product(p1, p6) *
               (scalar_product(p2, p4) - scalar_product(p3, p4) -
                scalar_product(p4, p5)) +
           (scalar_product(p1, p2) - scalar_product(p1, p3) -
            scalar_product(p1, p5)) *
               scalar_product(p4, p6) +
           scalar_product(p1, p4) *
               (scalar_product(p2, p6) - scalar_product(p3, p6) -
                scalar_product(p5, p6)))) /
             (-1 + scalar_product(p2, p2) - 2 * scalar_product(p2, p3) -
              2 * scalar_product(p2, p5) + scalar_product(p3, p3) +
              2 * scalar_product(p3, p5) + scalar_product(p5, p5)) +
         ((scalar_product(p2, p2) * scalar_product(p4, p5) -
           scalar_product(p2, p5) *
               (scalar_product(p4, p4) + 2 * scalar_product(p4, p5)) +
           scalar_product(p2, p4) *
               (2 * scalar_product(p2, p5) - 2 * scalar_product(p4, p5) -
                scalar_product(p5, p5))) *
          (scalar_product(p1, p6) *
               (scalar_product(p2, p3) - scalar_product(p3, p4) -
                scalar_product(p3, p5)) +
           (scalar_product(p1, p2) - scalar_product(p1, p4) -
            scalar_product(p1, p5)) *
               scalar_product(p3, p6) +
           scalar_product(p1, p3) *
               (scalar_product(p2, p6) - scalar_product(p4, p6) -
                scalar_product(p5, p6)))) /
             (-1 + scalar_product(p2, p2) - 2 * scalar_product(p2, p4) -
              2 * scalar_product(p2, p5) + scalar_product(p4, p4) +
              2 * scalar_product(p4, p5) + scalar_product(p5, p5)) +
         ((scalar_product(p3, p3) * scalar_product(p4, p5) +
           scalar_product(p3, p5) *
               (scalar_product(p4, p4) + 2 * scalar_product(p4, p5)) +
           scalar_product(p3, p4) *
               (2 * scalar_product(p3, p5) + 2 * scalar_product(p4, p5) +
                scalar_product(p5, p5))) *
          (scalar_product(p1, p6) *
               (scalar_product(p2, p3) + scalar_product(p2, p4) +
                scalar_product(p2, p5)) +
           scalar_product(p1, p3) * scalar_product(p2, p6) +
           scalar_product(p1, p4) * scalar_product(p2, p6) +
           scalar_product(p1, p5) * scalar_product(p2, p6) +
           scalar_product(p1, p2) * scalar_product(p3, p6) +
           scalar_product(p1, p2) * scalar_product(p4, p6) +
           scalar_product(p1, p2) * scalar_product(p5, p6))) /
             (-1 + scalar_product(p3, p3) + 2 * scalar_product(p3, p4) +
              2 * scalar_product(p3, p5) + scalar_product(p4, p4) +
              2 * scalar_product(p4, p5) + scalar_product(p5, p5)) +
         ((scalar_product(p2, p2) * scalar_product(p3, p4) -
           scalar_product(p2, p4) *
               (scalar_product(p3, p3) + 2 * scalar_product(p3, p4)) +
           scalar_product(p2, p3) *
               (2 * scalar_product(p2, p4) - 2 * scalar_product(p3, p4) -
                scalar_product(p4, p4))) *
          (scalar_product(p1, p6) *
               (scalar_product(p2, p5) - scalar_product(p3, p5) -
                scalar_product(p4, p5)) +
           scalar_product(p1, p5) *
               (scalar_product(p2, p6) - scalar_product(p3, p6) -
                scalar_product(p4, p6)) +
           (scalar_product(p1, p2) - scalar_product(p1, p3) -
            scalar_product(p1, p4)) *
               scalar_product(p5, p6))) /
             (-1 + scalar_product(p2, p2) - 2 * scalar_product(p2, p3) -
              2 * scalar_product(p2, p4) + scalar_product(p3, p3) +
              2 * scalar_product(p3, p4) + scalar_product(p4, p4)) +
         (((scalar_product(p1, p2) - scalar_product(p1, p3) -
            scalar_product(p1, p6)) *
               scalar_product(p4, p5) +
           scalar_product(p1, p5) *
               (scalar_product(p2, p4) - scalar_product(p3, p4) -
                scalar_product(p4, p6)) +
           scalar_product(p1, p4) *
               (scalar_product(p2, p5) - scalar_product(p3, p5) -
                scalar_product(p5, p6))) *
          (scalar_product(p2, p2) * scalar_product(p3, p6) -
           scalar_product(p2, p6) *
               (scalar_product(p3, p3) + 2 * scalar_product(p3, p6)) +
           scalar_product(p2, p3) *
               (2 * scalar_product(p2, p6) - 2 * scalar_product(p3, p6) -
                scalar_product(p6, p6)))) /
             (-1 + scalar_product(p2, p2) - 2 * scalar_product(p2, p3) -
              2 * scalar_product(p2, p6) + scalar_product(p3, p3) +
              2 * scalar_product(p3, p6) + scalar_product(p6, p6)) +
         (((scalar_product(p1, p2) - scalar_product(p1, p4) -
            scalar_product(p1, p6)) *
               scalar_product(p3, p5) +
           scalar_product(p1, p5) *
               (scalar_product(p2, p3) - scalar_product(p3, p4) -
                scalar_product(p3, p6)) +
           scalar_product(p1, p3) *
               (scalar_product(p2, p5) - scalar_product(p4, p5) -
                scalar_product(p5, p6))) *
          (scalar_product(p2, p2) * scalar_product(p4, p6) -
           scalar_product(p2, p6) *
               (scalar_product(p4, p4) + 2 * scalar_product(p4, p6)) +
           scalar_product(p2, p4) *
               (2 * scalar_product(p2, p6) - 2 * scalar_product(p4, p6) -
                scalar_product(p6, p6)))) /
             (-1 + scalar_product(p2, p2) - 2 * scalar_product(p2, p4) -
              2 * scalar_product(p2, p6) + scalar_product(p4, p4) +
              2 * scalar_product(p4, p6) + scalar_product(p6, p6)) +
         (((scalar_product(p1, p2) - scalar_product(p1, p5) -
            scalar_product(p1, p6)) *
               scalar_product(p3, p4) +
           scalar_product(p1, p4) *
               (scalar_product(p2, p3) - scalar_product(p3, p5) -
                scalar_product(p3, p6)) +
           scalar_product(p1, p3) *
               (scalar_product(p2, p4) - scalar_product(p4, p5) -
                scalar_product(p4, p6))) *
          (scalar_product(p2, p2) * scalar_product(p5, p6) -
           scalar_product(p2, p6) *
               (scalar_product(p5, p5) + 2 * scalar_product(p5, p6)) +
           scalar_product(p2, p5) *
               (2 * scalar_product(p2, p6) - 2 * scalar_product(p5, p6) -
                scalar_product(p6, p6)))) /
             (-1 + scalar_product(p2, p2) - 2 * scalar_product(p2, p5) -
              2 * scalar_product(p2, p6) + scalar_product(p5, p5) +
              2 * scalar_product(p5, p6) + scalar_product(p6, p6)) +
         ((scalar_product(p1, p3) * scalar_product(p2, p5) +
           scalar_product(p1, p4) * scalar_product(p2, p5) +
           scalar_product(p1, p6) * scalar_product(p2, p5) +
           scalar_product(p1, p5) *
               (scalar_product(p2, p3) + scalar_product(p2, p4) +
                scalar_product(p2, p6)) +
           scalar_product(p1, p2) * scalar_product(p3, p5) +
           scalar_product(p1, p2) * scalar_product(p4, p5) +
           scalar_product(p1, p2) * scalar_product(p5, p6)) *
          (scalar_product(p3, p3) * scalar_product(p4, p6) +
           scalar_product(p3, p6) *
               (scalar_product(p4, p4) + 2 * scalar_product(p4, p6)) +
           scalar_product(p3, p4) *
               (2 * scalar_product(p3, p6) + 2 * scalar_product(p4, p6) +
                scalar_product(p6, p6)))) /
             (-1 + scalar_product(p3, p3) + 2 * scalar_product(p3, p4) +
              2 * scalar_product(p3, p6) + scalar_product(p4, p4) +
              2 * scalar_product(p4, p6) + scalar_product(p6, p6)) +
         ((scalar_product(p1, p3) * scalar_product(p2, p4) +
           scalar_product(p1, p5) * scalar_product(p2, p4) +
           scalar_product(p1, p6) * scalar_product(p2, p4) +
           scalar_product(p1, p4) *
               (scalar_product(p2, p3) + scalar_product(p2, p5) +
                scalar_product(p2, p6)) +
           scalar_product(p1, p2) * scalar_product(p3, p4) +
           scalar_product(p1, p2) * scalar_product(p4, p5) +
           scalar_product(p1, p2) * scalar_product(p4, p6)) *
          (scalar_product(p3, p3) * scalar_product(p5, p6) +
           scalar_product(p3, p6) *
               (scalar_product(p5, p5) + 2 * scalar_product(p5, p6)) +
           scalar_product(p3, p5) *
               (2 * scalar_product(p3, p6) + 2 * scalar_product(p5, p6) +
                scalar_product(p6, p6)))) /
             (-1 + scalar_product(p3, p3) + 2 * scalar_product(p3, p5) +
              2 * scalar_product(p3, p6) + scalar_product(p5, p5) +
              2 * scalar_product(p5, p6) + scalar_product(p6, p6)) +
         ((scalar_product(p1, p4) * scalar_product(p2, p3) +
           scalar_product(p1, p5) * scalar_product(p2, p3) +
           scalar_product(p1, p6) * scalar_product(p2, p3) +
           scalar_product(p1, p3) * scalar_product(p2, p4) +
           scalar_product(p1, p3) * scalar_product(p2, p5) +
           scalar_product(p1, p3) * scalar_product(p2, p6) +
           scalar_product(p1, p2) * scalar_product(p3, p4) +
           scalar_product(p1, p2) * scalar_product(p3, p5) +
           scalar_product(p1, p2) * scalar_product(p3, p6)) *
          (scalar_product(p4, p4) * scalar_product(p5, p6) +
           scalar_product(p4, p6) *
               (scalar_product(p5, p5) + 2 * scalar_product(p5, p6)) +
           scalar_product(p4, p5) *
               (2 * scalar_product(p4, p6) + 2 * scalar_product(p5, p6) +
                scalar_product(p6, p6)))) /
             (-1 + scalar_product(p4, p4) + 2 * scalar_product(p4, p5) +
              2 * scalar_product(p4, p6) + scalar_product(p5, p5) +
              2 * scalar_product(p5, p6) + scalar_product(p6, p6));
}

//===========================================================================
//---- Amplitude for 2eta->4eta using only 6pt interactions -----------------
//===========================================================================
double amp_6pt(const FourMomentum &p1, const FourMomentum &p2,
               const FourMomentum &p3, const FourMomentum &p4,
               const FourMomentum &p5, const FourMomentum &p6) {
  return scalar_product(p1, p4) * scalar_product(p2, p6) *
             scalar_product(p3, p5) +
         scalar_product(p1, p4) * scalar_product(p2, p5) *
             scalar_product(p3, p6) +
         scalar_product(p1, p3) * scalar_product(p2, p6) *
             scalar_product(p4, p5) +
         scalar_product(p1, p2) * scalar_product(p3, p6) *
             scalar_product(p4, p5) +
         scalar_product(p1, p6) *
             (scalar_product(p2, p5) * scalar_product(p3, p4) +
              scalar_product(p2, p4) * scalar_product(p3, p5) +
              scalar_product(p2, p3) * scalar_product(p4, p5)) +
         scalar_product(p1, p3) * scalar_product(p2, p5) *
             scalar_product(p4, p6) +
         scalar_product(p1, p2) * scalar_product(p3, p5) *
             scalar_product(p4, p6) +
         scalar_product(p1, p5) *
             (scalar_product(p2, p6) * scalar_product(p3, p4) +
              scalar_product(p2, p4) * scalar_product(p3, p6) +
              scalar_product(p2, p3) * scalar_product(p4, p6)) +
         scalar_product(p1, p4) * scalar_product(p2, p3) *
             scalar_product(p5, p6) +
         scalar_product(p1, p3) * scalar_product(p2, p4) *
             scalar_product(p5, p6) +
         scalar_product(p1, p2) * scalar_product(p3, p4) *
             scalar_product(p5, p6);
}

std::tuple<double, double, double> scaled_cross_section(double z,
                                                        size_t nevents) {

  if (z <= 4.0) {
    return std::make_tuple(0.0, 0.0, 0.0);
  }

  // Incoming 4-momenta in the CM frame
  const double p = sqrt(z * z / 4.0 - 1.0);
  FourMomentum p1{z / 2.0, 0.0, 0.0, p};
  FourMomentum p2{z / 2.0, 0.0, 0.0, -p};

  auto msqrd44 = [&p1, &p2](const std::vector<FourMomentum> &fm) {
    const FourMomentum p3 = fm[0];
    const FourMomentum p4 = fm[1];
    const FourMomentum p5 = fm[2];
    const FourMomentum p6 = fm[3];

    const double amp = amp_4pt(p1, p2, p3, p4, p5, p6);
    return amp * amp;
  };
  auto msqrd66 = [&p1, &p2](const std::vector<FourMomentum> &fm) {
    const FourMomentum p3 = fm[0];
    const FourMomentum p4 = fm[1];
    const FourMomentum p5 = fm[2];
    const FourMomentum p6 = fm[3];

    const double amp = amp_6pt(p1, p2, p3, p4, p5, p6);
    return amp * amp;
  };
  auto msqrd46 = [&p1, &p2](const std::vector<FourMomentum> &fm) {
    const FourMomentum p3 = fm[0];
    const FourMomentum p4 = fm[1];
    const FourMomentum p5 = fm[2];
    const FourMomentum p6 = fm[3];

    const double amp4 = amp_4pt(p1, p2, p3, p4, p5, p6);
    const double amp6 = amp_6pt(p1, p2, p3, p4, p5, p6);

    return amp4 * amp6;
  };

  std::vector<double> isp_masses = {1.0, 1.0};
  std::vector<double> fsp_masses = {1.0, 1.0, 1.0, 1.0};

  auto res44 = Rambo(isp_masses, fsp_masses, z, msqrd44)
                   .compute_width_cross_section(nevents);
  auto res66 = Rambo(isp_masses, fsp_masses, z, msqrd66)
                   .compute_width_cross_section(nevents);
  auto res46 = Rambo(isp_masses, fsp_masses, z, msqrd46)
                   .compute_width_cross_section(nevents);
  // 24 is for a 4! since all FS particles are identical
  return std::make_tuple(res44.first / 24.0, res66.first / 24.0,
                         res46.first / 24.0);
}
