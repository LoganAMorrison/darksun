//
// Created by Logan Morrison on 2019-06-22.
//

#pragma once

#include "darksun/thermodynamic_particle.hpp"
#include "darksun/standard_model.hpp"
#include "darksun/compute_xi.hpp"
#include "darksun/odeint.hpp"
#include "darksun/stepper_sie.hpp"
#include <boost/math/differentiation/finite_difference.hpp>
#include <Eigen/Dense>

namespace darksun {

constexpr double PLANK_MASS = 1.220910e19;
constexpr double CRITICAL_ENERGY_DENSITY = 1.05375e-5;
constexpr double SM_ENTROPY_DENSITY_TODAY = 2891.2;
constexpr double T_CMB = 2.56215e-10;
constexpr double T_BBN = 0.0001; // 0.1 MeV in GeV

class DarkSUNModel : public ThermallyDecoupledModel {
protected:
private:
public:
    double lam;
    unsigned int N;
    double L1;
    double c;
    double a;
    double mu_eta;
    double mu_delta;
    bool has_dp;
    ThermodynamicParticle eta;
    ThermodynamicParticle delta;
    ThermodynamicParticle dark_photon;

    double xi_fo{-1};
    double Tsm_fo{-1};

    double xi_bbn{-1};
    double xi_cmb{-1};

    DarkSUNModel(double lam, unsigned int N, double L1, double c, double a,
                 double mu_eta, double mu_delta, double xi_inf_in, bool has_dp)
            : lam(lam), N(N), L1(L1), c(c), a(a), mu_eta(mu_eta), mu_delta(mu_delta), has_dp(has_dp),
              eta(ThermodynamicParticle{mu_eta * lam / sqrt(N), 1.0, 0}),
              delta(ThermodynamicParticle{mu_delta * lam * N, double(N + 1), N}),
              dark_photon(ThermodynamicParticle{0.0, has_dp ? 1.0 : 0.0, 2}) {

        xi_inf = xi_inf_in;
        hd_inf = 7.0 / 2.0 * N + 2 * N * N - 2.0 + dark_photon.g;
        sum_g = eta.g + dark_photon.g + (N % 2 == 0 ? 1.0 : 7.0 / 8.0) * delta.g;

        if (has_dp && dark_photon.mass < eta.mass) {
            gl = dark_photon.g;
            ml = dark_photon.mass;
        } else {
            gl = eta.g;
            ml = eta.mass;
        }
    }

    ~DarkSUNModel() = default;

    double dark_dof_entropy(double Td) const override {
        return eta.dof_entropy(Td) + delta.dof_entropy(Td) + dark_photon.dof_entropy(Td);
    }

    double dark_dof_energy(double Td) const {
        return eta.dof_energy(Td) + delta.dof_energy(Td) + dark_photon.dof_energy(Td);
    }

    double compute_xi(double Tsm) const {
        if (Tsm_fo != -1 || dark_photon.g > 0.0 || Tsm > Tsm_fo) {
            return compute_xi_const_tsm(Tsm);
        }
        if (Tsm_fo * xi_fo > eta.mass) return xi_fo;
        return xi_fo * Tsm / Tsm_fo;
    }

    double sqrt_gstar(double Tsm, double xi) const {
        double gd = dark_dof_energy(Tsm * xi);
        double gsm = sm_dof_energy(Tsm);
        return sm_sqrt_gstar(Tsm) * sqrt(gsm / (gsm + gd * xi * xi * xi * xi));
    }

    double thermal_cross_section_2eta_4eta(double x) const;
    double thermal_cross_section_2eta_2delta(double x) const;

    double cross_section_2eta_2eta() const;
    double cross_section_2delta_2delta() const;

    double delta_n_eff_cmb() const;
    double delta_n_eff_bbn() const;

    void operator()(double logx, Eigen::VectorXd &w, Eigen::VectorXd &dwdlogx) {
        double x = exp(logx);
        double Tsm = eta.mass / x;
        double xi = compute_xi(Tsm);
        double Td = xi * Tsm;
        double s = sm_entropy_density(Tsm);

        double neq = eta.neq(Td);
        double weq = log(neq / s);

        // Determine if the eta' has frozen out
        if (w(0) - weq > 0.4 && xi_fo < 0.0) {
            xi_fo = xi;
            Tsm_fo = Tsm;
        }
        // save xi at BBN
        if (Tsm < T_BBN && xi_bbn < 0.0) {
            xi_bbn = xi;
        }

        double sig_ee_eeee = thermal_cross_section_2eta_4eta(eta.mass / Td);
        double sig_eeee_ee = sig_ee_eeee / neq / neq;
        sig_eeee_ee = std::isnan(sig_eeee_ee) ? 0.0 : sig_eeee_ee;

        double sig_ee_dd = thermal_cross_section_2eta_2delta(eta.mass / Td);

        double pf_e = -sqrt(M_PI / 45) * PLANK_MASS * sqrt_gstar(Tsm, xi) * Tsm * s * s;
        double pf_d = sqrt(M_PI / 45) * PLANK_MASS * sqrt_gstar(Tsm, xi) * Tsm;

        // dW_e / dlogx
        dwdlogx(0) = pf_e * sig_eeee_ee * exp(w(0)) * (exp(2.0 * w(0)) - exp(2.0 * weq));

        // dY_d / dlogx
        dwdlogx(1) = pf_d * sig_ee_dd * exp(2.0 * w(0));
    }

    void jacobian(double logx, Eigen::VectorXd &w, Eigen::VectorXd &dfdx, Eigen::MatrixXd &dfdy) {
        using namespace boost::math::differentiation;

        double x = exp(logx);
        double Tsm = eta.mass / x;
        double xi = compute_xi(Tsm);
        double Td = xi * Tsm;
        double s = sm_entropy_density(Tsm);

        double neq = eta.neq(Td);
        double weq = log(neq / s);

        double sig_ee_eeee = thermal_cross_section_2eta_4eta(eta.mass / Td);
        double sig_eeee_ee = sig_ee_eeee / neq / neq;
        sig_eeee_ee = std::isnan(sig_eeee_ee) ? 0.0 : sig_eeee_ee;

        double sig_ee_dd = thermal_cross_section_2eta_2delta(eta.mass / Td);

        double pf_e = -sqrt(M_PI / 45) * PLANK_MASS * sqrt_gstar(Tsm, xi) * Tsm * s * s;
        double pf_d = sqrt(M_PI / 45) * PLANK_MASS * sqrt_gstar(Tsm, xi) * Tsm;

        dfdx(0) = finite_difference_derivative(
                [this, &w, dfdy](double logx) {
                    Eigen::VectorXd dw(w.size());
                    operator()(logx, w, dw);
                    return dw(0);
                }, logx);
        // df_d / dW_e
        dfdy(0, 0) = pf_e * sig_eeee_ee * exp(w(0)) * (3.0 * exp(2.0 * w(0)) - exp(2.0 * weq));

        // df_e / dY_d
        dfdy(0, 1) = 0.0;

        // df_d / dW_e
        dfdy(1, 0) = 2.0 * pf_d * sig_ee_dd * exp(2.0 * w(0));

        // df_d / dY_d
        dfdy(1, 1) = 0.0;
    }

    std::pair<double, double> compute_relic_density() {

        // reset parameters
        xi_fo = -1;
        Tsm_fo = -1;
        xi_bbn = -1;
        xi_cmb = -1;

        double Td_init = lam / 2.0;
        double xi_init = compute_xi_const_td(Td_init);
        double Tsm_init = Td_init / xi_init;

        double x_start = eta.mass / Tsm_init;
        double x_final = eta.mass / T_CMB;
        Eigen::VectorXd w_init(2);
        w_init(0) = log(eta.neq(Td_init) / sm_entropy_density(Tsm_init));

        double y_d_init = exp(-a * N) * delta.neq(Td_init) / sm_entropy_density(Tsm_init);
        w_init(1) = std::isnan(y_d_init) ? 0.0 : y_d_init;

        Output output(100);
        Odeint<StepperSie<DarkSUNModel>> solver(
                w_init, log(x_start), log(x_final), 1e-9, 1e-4, 1e-3, 0.0, output, *this);

        solver.integrate();

        double yinf_e = exp(output.ysave(0, output.count - 1));
        double yinf_d = output.ysave(1, output.count - 1);

        // save xi at CMB
        xi_cmb = compute_xi(T_CMB);

        return std::make_pair(
                yinf_e * eta.mass * SM_ENTROPY_DENSITY_TODAY / CRITICAL_ENERGY_DENSITY,
                yinf_d * delta.mass * SM_ENTROPY_DENSITY_TODAY / CRITICAL_ENERGY_DENSITY
        );
    }
};

}
