//
// Created by logan on 8/10/20.
//

#ifndef DARKSUN_SCANNER_HPP
#define DARKSUN_SCANNER_HPP

#include "darksun/darksun.hpp"
#include <fstream>
#include <functional>
#include <string>
#include <utility>

namespace darksun {

using ModelSetter = std::function<bool(size_t, DarkSunParameters &)>;

class Scanner {
public:
  const std::string file_name;
  ModelSetter set_model;

  Scanner(const std::string &t_file_name, ModelSetter t_set_model)
      : file_name(t_file_name), set_model(std::move(t_set_model)) {}

  void scan();

private:
  // Header for the output file
  const std::string header =
      "N,LAM,C,ADEL,LEC1,LEC2,MU_ETA,MU_DEL,XI_INF,XI_FO,TSM_FO,XI_CMB,XI_BBN,"
      "RD_ETA,RD_DEL,DNEFF_CMB,DNEFF_BBN,ETA_SI_PER_MASS,DEL_SI_PER_MASS";

  static void output_data(std::ofstream &ofile,
                          const DarkSunParameters &params);
};

void Scanner::scan() {
  size_t counter = 0;

  std::ofstream ofile;
  ofile.open(file_name);
  if (ofile.is_open()) {
    ofile << header;
  } else {
    throw std::runtime_error("Cannot open file: " + file_name);
  }

  while (true) {
    auto params = DarkSunParameters{0, 0};
    if (set_model(counter, params)) {
      break;
    }
    try {
      solve_boltzmann(1e-7, 1e-7, params);
    } catch (...) {
      params.xi_fo = NAN;
      params.tsm_fo = NAN;
      params.rd_eta = NAN;
      params.rd_del = NAN;
      params.xi_cmb = NAN;
      params.xi_bbn = NAN;
      params.dneff_cmb = NAN;
      params.dneff_bbn = NAN;
      params.eta_si_per_mass = NAN;
      params.del_si_per_mass = NAN;
    }
    output_data(ofile, params);
    counter++;
  }

  ofile.close();
}

void Scanner::output_data(std::ofstream &ofile,
                          const DarkSunParameters &params) {
  // Go to next line
  ofile << std::endl;
  // Model parameters
  ofile << params.n << ",";
  ofile << params.lam << ",";
  ofile << params.c << ",";
  ofile << params.adel << ",";
  ofile << params.lec1 << ",";
  ofile << params.lec2 << ",";
  ofile << params.mu_eta << ",";
  ofile << params.mu_del << ",";
  ofile << params.xi_inf << ",";
  // Derived parameters
  ofile << params.xi_fo << ",";
  ofile << params.tsm_fo << ",";
  ofile << params.xi_cmb << ",";
  ofile << params.xi_bbn << ",";
  ofile << params.rd_eta << ",";
  ofile << params.rd_del << ",";
  ofile << params.dneff_cmb << ",";
  ofile << params.dneff_bbn << ",";
  ofile << params.eta_si_per_mass << ",";
  ofile << params.del_si_per_mass;
}
} // namespace darksun

#endif // DARKSUN_SCANNER_HPP
