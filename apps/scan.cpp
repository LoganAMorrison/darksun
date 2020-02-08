//
// Created by Logan Morrison on 2/5/20.
//

#include "darksun/models.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <exception>
#include <boost/progress.hpp>

using namespace darksun;

void read_params(const std::string &file_name, double &L1, double &c, double &a,
                 double &mu_eta, double &mu_delta, double &xi_inf, bool &has_dp,
                 double &log_lam_min, double &log_lam_max, unsigned int &num_lams,
                 unsigned int &N_min, unsigned int &N_max, std::string &out_file);

int run_scan(double L1, double c, double a,
             double mu_eta, double mu_delta, double xi_inf, bool has_dp,
             double log_lam_min, double log_lam_max, unsigned int num_lams,
             unsigned int N_min, unsigned int N_max, const std::string &out_file);

int main(int argc, char *argv[]) {

    if (argc != 2) {
        std::cout << "Invalid number of arguments." << std::endl;
        std::cout << "Use: ./scan data.txt" << std::endl;
        return 1;
    } else {
        std::string file_name = argv[1];
        std::string out_file;
        double L1, c, a, mu_eta, mu_delta, xi_inf, log_lam_min, log_lam_max;
        unsigned int num_lams, N_min, N_max;
        bool has_dp = false;

        read_params(file_name, L1, c, a, mu_eta, mu_delta, xi_inf, has_dp,
                    log_lam_min, log_lam_max, num_lams, N_min, N_max,
                    out_file);

        return run_scan(L1, c, a, mu_eta, mu_delta, xi_inf, has_dp, log_lam_min,
                        log_lam_max, num_lams, N_min, N_max, out_file);

    }

}

void read_params(const std::string &file_name, double &L1, double &c, double &a,
                 double &mu_eta, double &mu_delta, double &xi_inf, bool &has_dp,
                 double &log_lam_min, double &log_lam_max, unsigned int &num_lams,
                 unsigned int &N_min, unsigned int &N_max, std::string &out_file) {

    int count = 0;
    std::string line;
    std::ifstream file(file_name);
    if (file.is_open()) {
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::vector<std::string> results(std::istream_iterator<std::string>{iss},
                                             std::istream_iterator<std::string>());

            if (results[0] == "L1") {
                L1 = std::stod(results[1]);
            } else if (results[0] == "c") {
                c = std::stod(results[1]);
            } else if (results[0] == "a") {
                a = std::stod(results[1]);
            } else if (results[0] == "mu_eta") {
                mu_eta = std::stod(results[1]);
            } else if (results[0] == "mu_delta") {
                mu_delta = std::stod(results[1]);
            } else if (results[0] == "xi_inf") {
                xi_inf = std::stod(results[1]);
            } else if (results[0] == "has_dp") {
                has_dp = bool(std::stoi(results[1]));
            } else if (results[0] == "log_lam_min") {
                log_lam_min = std::stod(results[1]);
            } else if (results[0] == "log_lam_max") {
                log_lam_max = std::stod(results[1]);
            } else if (results[0] == "num_lams") {
                num_lams = std::stoi(results[1]);
            } else if (results[0] == "N_max") {
                N_max = std::stoi(results[1]);
            } else if (results[0] == "N_min") {
                N_min = std::stoi(results[1]);
            } else if (results[0] == "out_file") {
                out_file = results[1];
            }
            count++;
        }
        file.close();
    }

    if (count != 13) {
        throw std::runtime_error("Invalid number of parameters in input file.");
    }

}

int run_scan(const double L1, const double c, const double a,
             const double mu_eta, const double mu_delta, const double xi_inf, const bool has_dp,
             const double log_lam_min, const double log_lam_max, const unsigned int num_lams,
             const unsigned int N_min, const unsigned int N_max, const std::string &out_file) {

    std::cout << "Running scan with the following parameters:" << std::endl;
    std::cout << "L1 = " << L1 << std::endl;
    std::cout << "c = " << c << std::endl;
    std::cout << "a = " << a << std::endl;
    std::cout << "mu_eta = " << mu_eta << std::endl;
    std::cout << "mu_delta = " << mu_delta << std::endl;
    std::cout << "xi_inf = " << xi_inf << std::endl;
    std::cout << "has_dp = " << has_dp << std::endl;
    std::cout << "log_lam_min = " << log_lam_min << std::endl;
    std::cout << "log_lam_max = " << log_lam_max << std::endl;
    std::cout << "num_lams = " << num_lams << std::endl;
    std::cout << "N_min = " << N_min << std::endl;
    std::cout << "N_max = " << N_max << std::endl;
    std::cout << "out_file = " << out_file << std::endl;


    std::vector<double> lams(num_lams);

    const double log_lam_step = (log_lam_max - log_lam_min) / double(num_lams - 1);
    for (size_t i = 0; i < num_lams; i++) {
        lams[i] = pow(10.0, log_lam_min + i * log_lam_step);
    }

    std::vector<unsigned int> Ns;
    Ns.reserve((N_max - N_min));
    for (unsigned int i = N_min; i <= N_max; i++) {
        Ns.push_back(i);
    }

    std::ofstream file;
    file.open(out_file);

    // Header
    file << "# L1 = " << L1 << ", c = " << c << ", a = " << a
         << ", mu_eta = " << mu_eta << ", mu_delta = " << mu_delta
         << ", xi_inf = " << xi_inf << ", has_dp = " << has_dp << std::endl;

    // Column names
    file << "LAMBDA, N, RD_ETA, RD_DELTA, (SIG_ETA / M_ETA), " <<
         "(SIG_DELTA / M_DELTA), DELTA_NEFF_CMB, DELTA_NEFF_BBN" <<
         std::endl;

    boost::progress_display progress(lams.size() * Ns.size());
    for (double lam : lams) {
        for (auto N: Ns) {
            try {
                DarkSUNModel model{lam, N, L1, c, a, mu_eta, mu_delta, xi_inf, has_dp};
                auto rds = model.compute_relic_density();
                double si_eta = model.cross_section_2eta_2eta() / model.eta.mass;
                double si_delta = model.cross_section_2delta_2delta() / model.delta.mass;
                double d_neff_cmb = model.delta_n_eff_cmb();
                double d_neff_bbn = model.delta_n_eff_bbn();

                file << lam << ", " << N << ", " << rds.first << ", " << rds.second << ", "
                     << si_eta << ", " << si_delta << ", "
                     << d_neff_cmb << ", " << d_neff_bbn << std::endl;
            } catch (...) {
                file << "nan" << ", " << "nan" << ", " << "nan" << ", " << "nan" << ", "
                     << "nan" << ", " << "nan" << ", "
                     << "nan" << ", " << "nan" << std::endl;
            }
            ++progress;
        }
    }

    file.close();

    return 0;
}