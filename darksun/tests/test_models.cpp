//
// Created by Logan Morrison on 2/5/20.
//

#include "darksun/models.hpp"
#include <iostream>

using namespace darksun;

int main() {

    constexpr double lam = 1e0;
    constexpr unsigned int N = 12;
    constexpr double L1 = 1.0;
    constexpr double c = 1.0;
    constexpr double a = 1.0;
    constexpr double mu_eta = 1.0;
    constexpr double mu_delta = 1.0;
    constexpr double xi_inf = 1e-1;
    constexpr bool has_dp = false;

    DarkSUNModel model{lam, N, L1, c, a, mu_eta, mu_delta, xi_inf, has_dp};

    auto rd = model.compute_relic_density();

    std::cout << "relic density eta = " << rd.first << std::endl;
    std::cout << "relic density delta = " << rd.second << std::endl;

    return 0;
}