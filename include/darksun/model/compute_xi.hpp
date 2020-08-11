#ifndef DARKSUN_MODEL_COMPUTE_XI_HPP
#define DARKSUN_MODEL_COMPUTE_XI_HPP

#include "darksun/model/parameters.hpp"
#include "darksun/model/thermal_functions.hpp"
#include <boost/math/tools/roots.hpp>
#include <gsl/gsl_sf_lambert.h>

namespace darksun {

auto xi_bounds_const_td(const double td, const DarkSunParameters &params)
    -> std::pair<double, double> {
  const double hd = dark_heff(td, params);
  const double cr_rhd = cbrt(dark_heff_inf(params) / hd);

  double lb = cr_rhd * cbrt(StandardModel::HEFF_0 / StandardModel::HEFF_INF) *
              params.xi_inf;
  double ub = cr_rhd * params.xi_inf;

  return std::make_pair(lb, ub);
}

auto xi_bounds_const_tsm(const double tsm, const DarkSunParameters &params)
    -> std::pair<double, double> {
  using boost::math::pow;
  const double xl = m_eta(params) / tsm;
  const double hsm = StandardModel::heff(tsm);
  const double hdinf = dark_heff_inf(params);
  const double sg = sum_g(params);

  const double lw_arg_num = pow<2>(45.0 * StandardModel::HEFF_INF * pow<3>(xl));
  const double lw_arg_den =
      pow<2>(4.0 * hdinf * hsm * params.xi_inf) * pow<7>(M_PI);

  const double ub = 2.0 * xl / gsl_sf_lambert_W0(lw_arg_num / lw_arg_den);
  const double lb =
      cbrt(hsm * hdinf / sg / StandardModel::HEFF_INF) * params.xi_inf;

  return std::make_pair(lb, ub);
}

double compute_xi_const_td(const double td, const DarkSunParameters &params) {
  using namespace boost::math;
  using namespace boost::math::tools;

  const double hd = dark_heff(td, params);
  const double c1 =
      dark_heff_inf(params) * pow<3>(params.xi_inf) / StandardModel::HEFF_INF;

  auto f = [hd, td, c1](double xi) -> double {
    return hd * pow<3>(xi) - StandardModel::heff(td / xi) * c1;
  };
  // Function specifying when to stop bisection algorithm
  auto tol = [](double min, double max) { return abs(max - min) <= 1e-8; };

  // Perform bisection algorith to find xi.
  const auto bounds = xi_bounds_const_td(td, params);
  const auto res = bisect(f, 0.8 * bounds.first, 1.2 * bounds.second, tol);
  // Return average of the bounding points
  return (res.second + res.first) / 2.0;
}

double compute_xi_const_tsm(const double tsm, const DarkSunParameters &params) {
  using namespace boost::math;
  using namespace boost::math::tools;

  if (params.tsm_fo >= 0.0) {
    if (tsm < params.tsm_fo) {
      // If tsm_fo has a value, that means the eta has frozen out. In this case,
      // there is no need to perform a root finding algorithm. We just
      // readshift.

      // Check if the eta is relativistic. If it is, its temperature redshifts
      // like the standard model temperature. Otherwise, it redshifts like
      // matter.
      if (params.tsm_fo * params.xi_fo > m_eta(params)) {
        return params.xi_fo;
      } else {
        return params.xi_fo * tsm / params.tsm_fo;
      }
    }
  }

  // The eta hasn't frozen out: assume that it is still in thermal
  // equillibrium with itself.
  const double hsm = StandardModel::heff(tsm);
  const double c1 = hsm * dark_heff_inf(params) * pow<3>(params.xi_inf) /
                    StandardModel::HEFF_INF;

  auto f = [c1, &params, tsm](double xi) -> double {
    return dark_heff(xi * tsm, params) * pow<3>(xi) - c1;
  };
  // Function specifying when to stop bisection algorithm
  auto tol = [](double min, double max) { return abs(max - min) <= 1e-8; };

  // Perform bisection algorith to find xi.
  const auto bounds = xi_bounds_const_tsm(tsm, params);
  const auto res = bisect(f, 0.8 * bounds.first, 1.2 * bounds.second, tol);
  // Return average of the bounding points
  return (res.second + res.first) / 2.0;
}

} // namespace darksun

#endif // DARKSUN_MODEL_COMPUTE_XI_HPP
