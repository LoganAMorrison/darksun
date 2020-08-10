#ifndef DARKSUN_MODEL_PARAMETERS_HPP
#define DARKSUN_MODEL_PARAMETERS_HPP

#include <cmath>
#include <cstdlib>
#include <gsl/gsl_spline.h>

namespace darksun {
namespace model {

class DarkSunParameters {
public:
  size_t n;             // N in SU(N)
  double lam;           // Confinement scale
  double c = 1.0;       // Suppression constant for 2eta->2del
  double adel = 1.0;    // Suppression constant for init delta abundance
  double lec1 = 0.1;    // Coefficient of 4pt eta interactions
  double lec2 = 1.0;    // Coefficient of 6pt eta interactions
  double mu_eta = 1.0;  // Coefficient of eta mass: me = mu * lam / sqrt(n)
  double mu_del = 1.0;  // Coefficient of del mass: md = mu * lam * n
  double xi_inf = 1e-2; // Ratio of dark to SM temperatures above EW scale

  double xi_fo = -1.0;
  double tsm_fo = -1.0;
  double xi_cmb = -1.0;
  double xi_bbn = -1.0;
  double rd_eta = -1.0;
  double rd_del = -1.0;
  double dneff_cmb = -1.0;
  double dneff_bbn = -1.0;
  double eta_si_per_mass = -1.0;
  double del_si_per_mass = -1.0;

  // Accelerators for use in interpolation function
  gsl_interp_accel *acc_cs44;
  gsl_interp_accel *acc_cs66;
  gsl_interp_accel *acc_cs46;

  double m_eta() const { return mu_eta * lam / sqrt(double(n)); }
  double m_del() const { return mu_del * lam * double(n); }
  double g_del() const { return double(n) + 1; }
  double dark_heff_inf() const { return 7.0 / 2.0 * n + 2.0 * (n * n - 1.0); }

  DarkSunParameters(size_t n, double lam) : n(n), lam(lam) {
    acc_cs44 = gsl_interp_accel_alloc();
    acc_cs66 = gsl_interp_accel_alloc();
    acc_cs46 = gsl_interp_accel_alloc();
  }
  ~DarkSunParameters() {
    gsl_interp_accel_free(acc_cs44);
    gsl_interp_accel_free(acc_cs66);
    gsl_interp_accel_free(acc_cs46);
  }

private:
};

double m_eta(const DarkSunParameters &params) {
  return params.mu_eta * params.lam / sqrt(double(params.n));
}
double m_del(const DarkSunParameters &params) {
  return params.mu_del * params.lam * double(params.n);
}
double g_del(const DarkSunParameters &params) { return double(params.n) + 1; }

double sum_g(const DarkSunParameters &params) { return 2.0 + double(params.n); }

} // namespace model
} // namespace darksun

#endif // DARKSUN_MODEL_PARAMETERS_HPP
