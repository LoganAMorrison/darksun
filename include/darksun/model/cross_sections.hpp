#ifndef DARKSUN_MODEL_CROSS_SECTIONS_HPP
#define DARKSUN_MODEL_CROSS_SECTIONS_HPP

#include "darksun/model/parameters.hpp"
#include "darksun/model/scaled_eta_cross_section.hpp"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <cmath>
#include <gsl/gsl_sf_bessel.h>

namespace darksun {
namespace model {

//===========================================================================
//---- Scaled cross-sections valid for all z = cme / meta -------------------
//===========================================================================

/**
 *  @breif Compute the cross-section for 2eta->4eta at zero temperature.
 *
 *  @param  cme  Center-of-mass energy.
 *  @param  model  Pointer to an array of doubles containing model parameters.
 *  @return Cross-section.
 *
 *  To compute the cross-section for 2eta->4eta, we use the scaled
 *  cross-sections (see 'DarkSun::scaled_cs_eta_24_xx') and add the model
 *  dependent prefactors.
 */
double cross_section_2eta_4eta(const double cme,
                               const DarkSunParameters &params) {

  using namespace boost::math;

  // Common prefactors. Constant is (256 pi^4 / 3)^2
  const double com = 6.9093374296577904e7 * pow<14>(params.mu_eta) /
                     pow<2>(params.lam) / pow<11>(params.n);
  // Coefficients of the A4*A4, A6*A6 and A4*A6 terms
  const double c44 = com * pow<4>(params.lec1) / 9.0;
  const double c66 = com * pow<2>(params.lec2) / 25.0;
  const double c46 = -2.0 * com * pow<2>(params.lec1) * params.lec2 / 15.0;
  // Scaled center-of-mass energy
  const double z = cme / params.m_eta();

  auto reff = ScaledEtaCrossSection::get_instance();

  return c44 * reff.scaled_cs_eta_44(z, params.acc_cs44) +
         c66 * reff.scaled_cs_eta_66(z, params.acc_cs66) +
         c46 * reff.scaled_cs_eta_46(z, params.acc_cs46);
}

/**
 *  @breif Compute the cross-section for 2eta->2eta at zero temperature.
 *
 *  @param  cme  Center-of-mass energy.
 *  @param  model  Pointer to an array of doubles containing model parameters.
 *  @return Cross-section.
 *
 *  The 2eta->2eta cross section is computed using the 4pt eta interaction
 *  from the chiral Lagrangian. It is assumed that the eta's are
 *  non-relativistic so that the energy of each eta is m_eta.
 */
double cross_section_2eta_2eta(const DarkSunParameters &params) {
  using namespace boost::math;
  return 62.012553360599640 * pow<6>(params.mu_eta) * pow<2>(params.lec1) /
         (pow<2>(params.lam) * pow<5>(params.n));
}

/**
 *  @breif Compute the cross-section for 2delta -> 2delta.
 *
 *  @param  cme  Center-of-mass energy.
 *  @param  model  Pointer to an array of doubles containing model parameters.
 *  @return Cross-section.
 *
 *  The 2delta->2delta cross section
 */
double cross_section_2del_2del(const DarkSunParameters &params) {
  using namespace boost::math;
  return 4.0 * pow<3>(M_PI) / pow<2>(params.lam);
}

//===========================================================================
//---- Thermally-averaged cross-sections ------------------------------------
//===========================================================================

double thermal_cross_section_2eta_4eta(const double x,
                                       const DarkSunParameters &params) {
  using boost::math::quadrature::gauss_kronrod;

  const double den = 2.0 * gsl_sf_bessel_Kn_scaled(2, x);
  const double pre = x / (den * den);
  const double meta = params.m_eta();

  auto f = [x, &params, meta](double z) -> double {
    const double z2 = z * z;
    const double sig = cross_section_2eta_4eta(z * meta, params);
    const double ker =
        z2 * (z2 - 4.0) * gsl_sf_bessel_K1_scaled(x * z) * exp(-x * (z - 2.0));
    return sig * ker;
  };

  const double integral = gauss_kronrod<double, 15>::integrate(
      f, 4.0, std::numeric_limits<double>::infinity(), 5, 1e-8);

  return pre * integral;
}

double thermal_cross_section_4eta_2eta(const double x,
                                       const DarkSunParameters &params) {
  using boost::math::pow;
  using boost::math::quadrature::gauss_kronrod;

  const double meta = params.m_eta();
  const double bes = gsl_sf_bessel_Kn_scaled(2, x);
  const double pre = pow<4>(M_PI) * pow<3>(x) / (pow<6>(meta) * pow<4>(bes));

  auto f = [x, &params, meta](double z) -> double {
    const double z2 = z * z;
    const double sig = cross_section_2eta_4eta(z * meta, params);
    const double ker =
        z2 * (z2 - 4.0) * gsl_sf_bessel_K1_scaled(x * z) * exp(-x * (z - 4.0));
    return sig * ker;
  };

  const double integral = gauss_kronrod<double, 15>::integrate(
      f, 4.0, std::numeric_limits<double>::infinity(), 5, 1e-8);

  return pre * integral;
}

double thermal_cross_section_2eta_2del(const double x,
                                       const DarkSunParameters &params) {
  using boost::math::pow;
  using boost::math::quadrature::gauss_kronrod;

  const double meta = params.m_eta();
  const double mdel = params.m_del();

  const double den = 2.0 * gsl_sf_bessel_Kn_scaled(2, x);
  const double pre = x / (den * den);
  const double zmin = 2.0 * mdel / meta;
  const double sig = exp(-2.0 * params.c * params.n) /
                     (64.0 * M_PI * pow<2>(params.n) * pow<2>(params.lam));

  auto f = [x](double z) -> double {
    const double z2 = z * z;
    return z2 * (z2 - 4.0) * gsl_sf_bessel_K1_scaled(x * z) *
           exp(-x * (z - 2.0));
  };

  const double integral = gauss_kronrod<double, 15>::integrate(
      f, zmin, std::numeric_limits<double>::infinity(), 15, 1e-8);

  return pre * sig * integral;
}

double thermal_cross_section_2del_2eta(const double xeta,
                                       const DarkSunParameters &params) {
  using boost::math::pow;
  using boost::math::quadrature::gauss_kronrod;

  const double meta = params.m_eta();
  const double mdel = params.m_del();
  const double xdel = mdel * xeta / meta;

  const double den = 2.0 * gsl_sf_bessel_Kn_scaled(2, xdel);
  const double pre = xdel / (den * den);
  const double sig = exp(-2.0 * params.c * params.n) /
                     (64.0 * M_PI * pow<2>(params.n) * pow<2>(params.lam));

  auto f = [xdel](double z) -> double {
    const double z2 = z * z;
    return z2 * (z2 - 4.0) * gsl_sf_bessel_K1_scaled(xdel * z) *
           exp(-xdel * (z - 2.0));
  };

  const double integral = gauss_kronrod<double, 15>::integrate(
      f, 2.0, std::numeric_limits<double>::infinity(), 15, 1e-8);

  return pre * sig * integral;
}

} // namespace model
} // namespace darksun

#endif // DARKSUN_MODEL_CROSS_SECTIONS_HPP
