#ifndef DARKSUN_MODEL_PARAMETERS_HPP
#define DARKSUN_MODEL_PARAMETERS_HPP

#include <cstdlib>
#include <gsl/gsl_spline.h>

namespace darksun {
namespace model {

class DarkSun {
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

  DarkSun(size_t n, double lam) : n(n), lam(lam) {}

private:
};

} // namespace model
} // namespace darksun

#endif // DARKSUN_MODEL_PARAMETERS_HPP
