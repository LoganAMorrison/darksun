#ifndef DARKSUN_MODEL_DNEFF_HPP
#define DARKSUN_MODEL_DNEFF_HPP

#include "darksun/constants.hpp"
#include "darksun/model/parameters.hpp"
#include "darksun/model/thermal_functions.hpp"
#include <boost/math/special_functions/pow.hpp>

namespace darksun {
namespace model {

double compute_dneff_cmb(const DarkSunParameters &params) {
  using namespace boost::math;
  const double xi = params.xi_cmb;
  if (xi >= 0.0) {

    return 4.0 / 7.0 * pow(11.0 / 4.0, 4.0 / 3.0) *
           dark_geff(T_CMB * xi, params) * pow<4>(xi);
  }
  return NAN;
}

double compute_dneff_bbn(const DarkSunParameters &params) {
  using namespace boost::math;
  const double xi = params.xi_bbn;
  if (xi >= 0.0) {
    return 4.0 / 7.0 * dark_geff(T_BBN * xi, params) * pow<4>(xi);
  }
  return NAN;
}
} // namespace model
} // namespace darksun

#endif // DARKSUN_MODEL_DNEFF_HPP
