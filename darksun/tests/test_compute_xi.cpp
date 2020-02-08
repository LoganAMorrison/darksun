//
// Created by Logan Morrison on 2/4/20.
//

#include "darksun/compute_xi.hpp"
#include "darksun/thermodynamic_particle.hpp"
#include "darksun/standard_model.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <cmath>

using namespace darksun;

class TestThermallyDecoupledModel : public ThermallyDecoupledModel {
public:
    ThermodynamicParticle eta;
    ThermodynamicParticle delta;
    ThermodynamicParticle dark_photon;

    TestThermallyDecoupledModel(double lam, unsigned int N, double xi_inf, bool has_dp)
            : eta(ThermodynamicParticle{lam / sqrt(N), 1.0, 0}),
              delta(ThermodynamicParticle{lam * N, double(N + 1), N}),
              dark_photon(ThermodynamicParticle{0.0, (has_dp ? 2.0 : 0.0), 1}) {
        this->xi_inf = xi_inf;
        // quarks + gluons + dark photon
        hd_inf = (4.0 * N * 7.0 / 8.0) + (2.0 * (N * N - 1.0)) + dark_photon.g;
        sum_g = eta.g + dark_photon.g + delta.g * (delta.spin2 % 2 == 0 ? 1.0 : 7.0 / 8.0);
        gl = (has_dp ? dark_photon.g : eta.g);
        ml = (has_dp ? 0.0 : eta.mass);
    }

    double dark_dof_entropy(double Td) const override {
        return eta.dof_entropy(Td) + delta.dof_entropy(Td) + dark_photon.dof_entropy(Td);
    }
};


TEST(ComputeXiTest, TestConstantTDark) {

    constexpr double log_Td_min = -3.0;
    constexpr double log_Td_max = 3.0;
    constexpr int num_Tds = 1000;
    const double log_Td_step = (log_Td_max - log_Td_min) / double(num_Tds - 1);

    for (int i = 0; i < num_Tds; i++) {
        double Td = pow(10.0, log_Td_min + i * log_Td_step);
        auto model = TestThermallyDecoupledModel{1.0, 50, 1.0, true};

        double xi = model.compute_xi_const_td(Td);
        double Tsm = Td / xi;

        double hsm = sm_dof_entropy(Tsm);
        double hd = model.dark_dof_entropy(Td);

        double residual = pow(xi, 3) - (hsm * model.hd_inf / hd / _hsm_inf) * pow(model.xi_inf, 3);

        ASSERT_LE(fabs(residual), 1e-4);
    }

    for (int i = 0; i < num_Tds; i++) {
        double Td = pow(10.0, log_Td_min + i * log_Td_step);
        auto model = TestThermallyDecoupledModel{1e0, 50, 1e0, false};

        double xi = model.compute_xi_const_td(Td);
        double Tsm = Td / xi;

        double hsm = sm_dof_entropy(Tsm);
        double hd = model.dark_dof_entropy(Td);

        double rhs = (hsm * model.hd_inf / hd / _hsm_inf) * pow(model.xi_inf, 3);
        double lhs = pow(xi, 3);

        double frac_diff = fabs(rhs - lhs) / rhs;

        ASSERT_LE(frac_diff, 1e-4);
    }
}

TEST(ComputeXiTest, TestConstantTDarkSmallXi) {

    constexpr double log_Td_min = -3.0;
    constexpr double log_Td_max = 3.0;
    constexpr int num_Tds = 1000;
    const double log_Td_step = (log_Td_max - log_Td_min) / double(num_Tds - 1);

    for (int i = 0; i < num_Tds; i++) {
        double Td = pow(10.0, log_Td_min + i * log_Td_step);
        auto model = TestThermallyDecoupledModel{1.0, 50, 1e-3, true};

        double xi = model.compute_xi_const_td(Td);
        double Tsm = Td / xi;

        double hsm = sm_dof_entropy(Tsm);
        double hd = model.dark_dof_entropy(Td);

        double residual = pow(xi, 3) - (hsm * model.hd_inf / hd / _hsm_inf) * pow(model.xi_inf, 3);

        ASSERT_LE(fabs(residual), 1e-4);
    }

    for (int i = 0; i < num_Tds; i++) {
        double Td = pow(10.0, log_Td_min + i * log_Td_step);
        auto model = TestThermallyDecoupledModel{1e0, 50, 1e-3, false};

        double xi = model.compute_xi_const_td(Td);
        double Tsm = Td / xi;

        double hsm = sm_dof_entropy(Tsm);
        double hd = model.dark_dof_entropy(Td);

        double rhs = (hsm * model.hd_inf / hd / _hsm_inf) * pow(model.xi_inf, 3);
        double lhs = pow(xi, 3);

        double frac_diff = fabs(rhs - lhs) / rhs;

        ASSERT_LE(frac_diff, 1e-4);
    }
}

TEST(ComputeXiTest, TestConstantTSM) {

    constexpr double log_Tsm_min = -3.0;
    constexpr double log_Tsm_max = 3.0;
    constexpr int num_Tsms = 1000;
    const double log_Tsm_step = (log_Tsm_max - log_Tsm_min) / double(num_Tsms - 1);

    for (int i = 0; i < num_Tsms; i++) {
        double Tsm = pow(10.0, log_Tsm_min + i * log_Tsm_step);
        auto model = TestThermallyDecoupledModel{1.0, 50, 1.0, true};

        double xi = model.compute_xi_const_tsm(Tsm);
        double Td = Tsm * xi;

        double hsm = sm_dof_entropy(Tsm);
        double hd = model.dark_dof_entropy(Td);

        double rhs = (hsm * model.hd_inf / hd / _hsm_inf) * pow(model.xi_inf, 3);
        double lhs = pow(xi, 3);

        double frac_diff = fabs(rhs - lhs) / rhs;

        ASSERT_LE(frac_diff, 1e-4);
    }

    for (int i = 0; i < num_Tsms; i++) {
        double Tsm = pow(10.0, log_Tsm_min + i * log_Tsm_step);
        auto model = TestThermallyDecoupledModel{1e0, 50, 1e0, false};

        double xi = model.compute_xi_const_tsm(Tsm);
        double Td = Tsm * xi;

        double hsm = sm_dof_entropy(Tsm);
        double hd = model.dark_dof_entropy(Td);

        double rhs = (hsm * model.hd_inf / _hsm_inf) * pow(model.xi_inf, 3);
        double lhs = hd * pow(xi, 3);

        double frac_diff = fabs(rhs - lhs) / rhs;

        ASSERT_LE(frac_diff, 1e-4);
    }
}

TEST(ComputeXiTest, TestConstantTSMSmallXi) {

    constexpr double log_Tsm_min = -3.0;
    constexpr double log_Tsm_max = 3.0;
    constexpr int num_Tsms = 1000;
    const double log_Tsm_step = (log_Tsm_max - log_Tsm_min) / double(num_Tsms - 1);

    for (int i = 0; i < num_Tsms; i++) {
        double Tsm = pow(10.0, log_Tsm_min + i * log_Tsm_step);
        auto model = TestThermallyDecoupledModel{1.0, 50, 1e-3, true};

        double xi = model.compute_xi_const_tsm(Tsm);
        double Td = Tsm * xi;

        double hsm = sm_dof_entropy(Tsm);
        double hd = model.dark_dof_entropy(Td);

        double rhs = (hsm * model.hd_inf / hd / _hsm_inf) * pow(model.xi_inf, 3);
        double lhs = pow(xi, 3);

        double frac_diff = fabs(rhs - lhs) / rhs;

        ASSERT_LE(frac_diff, 1e-4);
    }

    for (int i = 0; i < num_Tsms; i++) {
        double Tsm = pow(10.0, log_Tsm_min + i * log_Tsm_step);
        auto model = TestThermallyDecoupledModel{1e0, 50, 1e-3, false};

        double xi = model.compute_xi_const_tsm(Tsm);
        double Td = Tsm * xi;

        double hsm = sm_dof_entropy(Tsm);
        double hd = model.dark_dof_entropy(Td);

        double rhs = (hsm * model.hd_inf / _hsm_inf) * pow(model.xi_inf, 3);
        double lhs = hd * pow(xi, 3);

        double frac_diff = fabs(rhs - lhs) / rhs;

        ASSERT_LE(frac_diff, 1e-4);
    }
}


int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
