//
// Created by Logan Morrison on 2/7/20.
//

#include "darksun/compute_xi.hpp"
#include "darksun/standard_model.hpp"

namespace darksun {

/**
 * Lower bound on xi = Tsm / Td assuming Td is a constant.
 * @param Td Dark sector temperature.
 * @return Lower bound on xi.
 */
double ThermallyDecoupledModel::xi_lower_bound_const_td(double Td) const {
    double hd = dark_dof_entropy(Td);
    return cbrt(hd_inf / hd * _hsm0 / _hsm_inf) * xi_inf;
}

/**
 * Upper bound on xi = Tsm / Td assuming Td is a constant.
 * @param Td Dark sector temperature.
 * @return Upper bound on xi.
 */
double ThermallyDecoupledModel::xi_upper_bound_const_td(double Td) const {
    double hd = dark_dof_entropy(Td);
    return cbrt(hd_inf / hd) * xi_inf;
}

/**
 * Lower bound on xi = Tsm / Td assuming Tsm is a constant.
 * @param Tsm Temperature of the SM bath.
 * @return Lower bound.
 */
double ThermallyDecoupledModel::xi_lower_bound_const_tsm(double Tsm) const {
    return cbrt(sm_dof_entropy(Tsm) * hd_inf / sum_g / _hsm_inf) * xi_inf;
}

/**
 * Upper bound on xi = Tsm / Td assuming Tsm is a constant
 * @param Tsm Temperature of the SM bath.
 * @return Lower bound.
 */
double ThermallyDecoupledModel::xi_upper_bound_const_tsm(double Tsm) const {
    if (ml <= 0.0)
        return cbrt(sm_dof_entropy(Tsm) * hd_inf / gl / _hsm_inf) * xi_inf;
    else {
        using namespace boost::math;

        double xl = ml / Tsm;
        double hsm = sm_dof_entropy(Tsm);

        double lw_arg_num = pow(45.0 * gl * _hsm_inf * xl * xl * xl, 2);
        double lw_arg_den = pow(4.0 * hd_inf * hsm * xi_inf, 2) * pow(M_PI, 7);
        return 2.0 * xl / lambert_w0(lw_arg_num / lw_arg_den);
    }
}

/**
 * Compute the value of xi = Td / Tsm assuming Td is a constant.
 * @param Td Temperature of the dark bath.
 * @return Ratio of dark to SM temperatures: xi = Td / Tsm.
 */
double ThermallyDecoupledModel::compute_xi_const_td(double Td) const {
    using namespace boost::math::tools;

    double hd = dark_dof_entropy(Td);

    double lb = xi_lower_bound_const_td(Td);
    double ub = xi_upper_bound_const_td(Td);

    double c1 = hd_inf * pow(xi_inf, 3) / _hsm_inf;
    auto f = [hd, Td, c1](double xi) {
        return hd * pow(xi, 3) - sm_dof_entropy(Td / xi) * c1;
    };
    auto tol = [](double min, double max) { return abs(max - min) <= 1e-8; };

    std::pair<double, double> res = bisect(f, 0.8 * lb, 1.2 * ub, tol);
    return (res.second + res.first) / 2.0;
}

/**
 * Compute the value of xi = Td / Tsm assuming Tsm is a constant.
 * @param Tsm Temperature of the SM bath.
 * @return Ratio of dark to SM temperatures: xi = Td / Tsm.
 */
double ThermallyDecoupledModel::compute_xi_const_tsm(double Tsm) const {
    using namespace boost::math::tools;
    double ub = xi_upper_bound_const_tsm(Tsm);
    double lb = xi_lower_bound_const_tsm(Tsm);
    double hsm = sm_dof_entropy(Tsm);

    double c1 = hsm * hd_inf * pow(xi_inf, 3) / _hsm_inf;
    auto f = [this, Tsm, c1](double xi) {
        return dark_dof_entropy(xi * Tsm) * pow(xi, 3) - c1;
    };
    auto tol = [](double min, double max) { return abs(max - min) <= 1e-8; };

    std::pair<double, double> res = bisect(f, 0.8 * lb, 1.2 * ub, tol);
    return (res.second + res.first) / 2.0;
}

}