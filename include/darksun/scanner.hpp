//
// Created by logan on 8/10/20.
//

#ifndef DARKSUN_SCANNER_HPP
#define DARKSUN_SCANNER_HPP

#include "darksun/darksun.hpp"
#include <fstream>
#include <functional>
#include <mutex>
#include <string>
#include <thread>
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
  // Mutex for outputting data to file
  std::mutex outmutex;
  // Mutex for getting iter
  std::mutex itermutex;
  // Stream object for outputting data
  std::ofstream ofile;

  size_t iter = 0;
  size_t get_iter();
  void thread_scan();

  // Spawner for threads
  std::thread spawn_thread_scan() {
    return std::thread([this] { this->thread_scan(); });
  }

  void output_data(std::ofstream &ofile, const DarkSunParameters &params);
};

size_t Scanner::get_iter() {
  // Lock function, get current iter, update iter and return value
  std::lock_guard<std::mutex> lock(itermutex);
  size_t it = iter;
  iter++;
  return it;
}

void Scanner::thread_scan() {
  while (true) {
    auto params = DarkSunParameters{0, 0};
    size_t it = get_iter();
    if (set_model(it, params)) {
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
  }
}

void Scanner::scan() {
  // Reset the iter
  iter = 0;

  // Open and make sure output file exists. If it does, write the header.
  ofile.open(file_name);
  if (ofile.is_open()) {
    ofile << header;
  } else {
    throw std::runtime_error("Cannot open file: " + file_name);
  }

  // Function for threads to run
  // auto f = [this]() { thread_scan(); };

  // Get the number of cpus
  const auto cpu_count = std::thread::hardware_concurrency();
  // Create threads
  std::vector<std::thread> threads(cpu_count);
  // Launch the threads
  for (auto &thread : threads) {
    thread = spawn_thread_scan();
  }
  // Wait for threads to finish
  for (auto &thread : threads) {
    thread.join();
  }

  ofile.close();
}

void Scanner::output_data(std::ofstream &ofile,
                          const DarkSunParameters &params) {
  // Lock this function to avoid data races
  const std::lock_guard<std::mutex> lock(outmutex);
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
