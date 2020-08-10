#ifndef DARKSUN_MODEL_THERMAL_FUNCTIONS_HPP
#define DARKSUN_MODEL_THERMAL_FUNCTIONS_HPP

#include "darksun/model/parameters.hpp"
#include "darksun/standard_model.hpp"
#include <boost/math/special_functions/pow.hpp>
#include <gsl/gsl_sf_bessel.h>

namespace darksun {
namespace model {

//===========================================================================
//---- Functions for computing equilibrium number densities -----------------
//===========================================================================

double neq_eta(const double td, const DarkSunParameters &params) {
  using boost::math::pow;

  const double x = params.m_eta() / td;
  const double g = 1.0;
  const double fac = g * pow<2>(x) * pow<3>(td) / (2.0 * pow<2>(M_PI));

  double bessum = 0.0;
  for (int k = 0; k < 5; k++) {
    bessum += gsl_sf_bessel_Kn(2, (1 + k) * x) / (1 + k);
  }
  return fac * bessum;
}

double neq_del(const double td, const DarkSunParameters &params) {
  using boost::math::pow;

  const double x = params.m_del() / td;
  const double g = params.g_del();
  const double eta = params.n % 2 == 0 ? 1.0 : -1.0;
  const double fac = g * pow<2>(x) * pow<3>(td) / (2.0 * pow<2>(M_PI));

  double bessum = 0.0;
  for (int k = 0; k < 5; k++) {
    bessum += std::pow(eta, k) * gsl_sf_bessel_Kn(2, (1 + k) * x) / (1 + k);
  }
  return fac * bessum;
}

double yeq_eta(const double tsm, const double xi,
               const DarkSunParameters &params) {
  using boost::math::pow;

  const double x = params.m_eta() / (tsm * xi);
  const double g = 1.0;
  const double fac =
      45.0 * g * x * x / (4.0 * pow<4>(M_PI) * StandardModel::heff(tsm));

  double bessum = 0.0;
  for (int k = 0; k < 5; k++) {
    bessum += gsl_sf_bessel_Kn(2, (1 + k) * x) / (1 + k);
  }
  return fac * bessum;
}

double yeq_del(const double tsm, const double xi,
               const DarkSunParameters &params) {
  using boost::math::pow;
  const double x = params.m_del() / (tsm * xi);
  const double g = params.g_del();
  const double eta = params.n % 2 == 0 ? 1.0 : -1.0;
  const double fac =
      45.0 * g * x * x / (4.0 * pow<4>(M_PI) * StandardModel::heff(tsm));

  double bessum = 0.0;
  for (int k = 0; k < 5; k++) {
    bessum += std::pow(eta, k) * gsl_sf_bessel_Kn(2, (1 + k) * x) / (1 + k);
  }
  return fac * bessum;
}

double weq_eta(const double tsm, const double xi,
               const DarkSunParameters &params) {
  using boost::math::pow;
  const double x = params.m_eta() / (tsm * xi);
  const double g = 1.0;
  const double fac =
      45.0 * g * x * x / (4.0 * pow<4>(M_PI) * StandardModel::heff(tsm));

  double bessum = 0.0;
  for (int k = 0; k < 5; k++) {
    bessum += exp(-k * x) * gsl_sf_bessel_Kn_scaled(2, (1 + k) * x) / (1 + k);
  }
  return -x + log(fac * bessum);
}

double weq_del(const double tsm, const double xi,
               const DarkSunParameters &params) {
  using boost::math::pow;
  const double x = params.m_del() / (tsm * xi);
  const double g = params.g_del();
  const double eta = params.n % 2 == 0 ? 1.0 : -1.0;
  const double fac =
      45.0 * g * x * x / (4.0 * pow<4>(M_PI) * StandardModel::heff(tsm));

  double bessum = 0.0;
  for (int k = 0; k < 5; k++) {
    bessum += exp(-k * x) * std::pow(eta, k) *
              gsl_sf_bessel_Kn_scaled(2, (1 + k) * x) / (1 + k);
  }
  return -x + log(fac * bessum);
}

//===========================================================================
//---- Misc. thermal functions ----------------------------------------------
//===========================================================================

double dark_heff_inf(const DarkSunParameters &params) {
  return 7.0 / 2.0 * params.n + 2.0 * (params.n * params.n - 1.0);
}

double dark_heff(const double td, const DarkSunParameters &params) {
  using boost::math::pow;
  const double xe = params.m_eta() / td;
  const double xd = params.m_del() / td;

  const double ge = 1.0;
  const double gd = params.g_del();

  const double etae = 1.0;
  const double etad = params.n % 2 == 0 ? 1.0 : -1.0;

  const double pre = 45.0 / (4.0 * pow<4>(M_PI));
  const double pree = pre * ge * pow<3>(xe);
  const double pred = pre * gd * pow<3>(xd);

  double bessume = 0.0;
  double bessumd = 0.0;
  for (int k = 0; k < 5; k++) {
    bessume += std::pow(etae, k) / (1 + k) * gsl_sf_bessel_Kn(3, (1 + k) * xe);
    bessumd += std::pow(etad, k) / (1 + k) * gsl_sf_bessel_Kn(3, (1 + k) * xd);
  }

  return pree * bessume + pred * bessumd;
}

double dark_geff(const double td, const DarkSunParameters &params) {
  using boost::math::pow;
  const double xe = params.m_eta() / td;
  const double xd = params.m_del() / td;

  const double ge = 1.0;
  const double gd = params.g_del();

  const double etae = 1.0;
  const double etad = params.n % 2 == 0 ? 1.0 : -1.0;

  const double pre = 30.0 / (2.0 * pow<4>(M_PI));
  const double pree = pre * ge * pow<2>(xe);
  const double pred = pre * gd * pow<2>(xd);

  double bessume = 0.0;
  double bessumd = 0.0;
  for (int k = 0; k < 5; k++) {
    bessume += std::pow(etae, k) / pow<2>(1 + k) *
               ((1 + k) * xe * gsl_sf_bessel_Kn(1, (1 + k) * xe) +
                3.0 * gsl_sf_bessel_Kn(2, (1 + k) * xe));
    bessumd += std::pow(etad, k) / pow<2>(1 + k) *
               ((1 + k) * xd * gsl_sf_bessel_Kn(1, (1 + k) * xd) +
                3.0 * gsl_sf_bessel_Kn(2, (1 + k) * xd));
  }

  return pree * bessume + pred * bessumd;
}

double sqrt_gstar(const double tsm, const double xi,
                  const DarkSunParameters &params) {
  using boost::math::pow;
  double gd = dark_geff(tsm * xi, params);
  double gsm = StandardModel::geff(tsm);
  return StandardModel::sqrt_gstar(tsm) * sqrt(gsm / (gsm + gd * pow<4>(xi)));
}

} // namespace model
} // namespace darksun

#endif // DARKSUN_MODEL_THERMAL_FUNCTIONS_HPP
