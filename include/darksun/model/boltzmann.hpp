#ifndef DARKSUN_MODEL_BOLTZMANN_HPP
#define DARKSUN_MODEL_BOLTZMANN_HPP

#include "darksun/model/compute_xi.hpp"
#include "darksun/model/parameters.hpp"
#include "darksun/model/thermal_functions.hpp"
#include <gsl/gsl_errno.h>

namespace darksun {
namespace model {

void solve_boltzmann(double reltol, double abstol, DarkSunParameters &params) {
  gsl_set_error_handler_off();

  // Initial conditions
  double meta = params.m_eta();

  // Integration interval
  double td = params.lam / 2.0;                // Start Td at confinement
  double xi = compute_xi_const_td(td, params); // Starting value of xi
  double tsm = td / xi;                        // Initial SM temperature
  double start = log(meta / tsm);
  double final = log(meta / T_CMB);

  // Initial conditions : y[0] = log(Y_eta), y[1] = Y_del
  constexpr int ndim = 2;
  double y[ndim];
  y[0] = weq_eta(tsm, xi, params);
  y[1] = exp(-params.adel * params.n) * yeq_del(tsm, xi, params);

  //==================================================================
  //---- Set RADAU scalar parameters ---------------------------------
  //==================================================================
  int nd = ndim;        // dimension of the system
  int ns = 9;           // maximum number of allowed stages
  int ijac = 1;         // analytic jacobian
  int mljac = ndim;     // full Jacaobian
  int mujac = 0;        // full Jacaobian
  int imas = 0;         // no mass matrix
  int mlmas = 0;        // no mass matrix
  int mumas = 0;        // no mass matrix
  int idid{};           // Flag for radau return code
  double rtol = reltol; // Relative tolerance for radau
  double atol = abstol; // Absolute tolerance for radau
  int itol = 0;         // Flag telling radau we are using scalar tols
  double h = 1.0e-6;    // Initial step size
  int iout = 1;         // Have radau call `solout`
  int ipar = 0;         // Counter for specifying solution index
  int lwork = (ns + 1) * nd * nd + (3 * ns + 3) * nd + 20;
  int liwork = (2 + (ns - 1) / 2) * nd + 20;

  //==================================================================
  //---- Set RADAU array parameters ----------------------------------
  //==================================================================
  std::vector<double> work(lwork); // Workspace of doubles needed for radau
  std::vector<int> iwork(liwork);  // Workspace of ints needed for radau
  for (int i = 0; i < 20; i++) {   // Set all radau params to defaults
    iwork[i] = 0;
    work[i] = 0.0;
  }

  // Set the spacing between solutions
  set_dx(model, (final - start) / double(DarkSun::NUM_SOLS - 1));

  //==================================================================
  //---- Define lambdas which capture the model ----------------------
  //==================================================================
  auto boltz = [model, &ipar](int *n, double *logx, double *y, double *dy) {
    DarkSun::boltzmann(n, logx, y, dy, model, &ipar);
  };
  auto jac = [model, &ipar](int *n, double *logx, double *y, double *dfy,
                            int *ldfy) {
    DarkSun::boltzmann_jac(n, logx, y, dfy, ldfy, model, &ipar);
  };
  auto mas = [](int *, double *, int *) {};
  auto solo = [model, &ipar](int *nr, double *xold, double *x, double *y,
                             double *cont, int *lrc, int *n, int *irtrn,
                             const RadauWeight &w) {
    DarkSun::solout(nr, xold, x, y, cont, lrc, n, model, &ipar, irtrn, w);
  };

  //==================================================================
  //---- Solve the Boltzmann equations using RADAU -------------------
  //==================================================================
  radau(&nd, boltz, &start, y, &final, &h, &rtol, &atol, &itol, jac, &ijac,
        &mljac, &mujac, mas, &imas, &mlmas, &mumas, solo, &iout, work.data(),
        &lwork, iwork.data(), &liwork, &idid);

  // Save last state
  set_sol_boltz(DarkSun::NUM_SOLS - 1, final, y[0], y[1], model);

  //==================================================================
  //---- Record/Calculate outputs ------------------------------------
  //==================================================================
  // If 'idid' > 0, that means that RADAU was successful. Otherwise,
  // there was an error
  if (idid > 0) {
    set_rd_eta(model, m_eta(model) * exp(y[0]) * S_TODAY / RHO_CRIT);
    set_rd_del(model, m_del(model) * y[1] * S_TODAY / RHO_CRIT);
    set_xi_cmb(model, compute_xi_const_tsm(T_CMB, model));
    set_xi_bbn(model, compute_xi_const_tsm(T_BBN, model));
    set_dneff_bbn(model, compute_dneff_bbn(model));
    set_dneff_cmb(model, compute_dneff_cmb(model));
    set_eta_si_per_mass(model, cross_section_2eta_2eta(model) / m_eta(model));
    set_del_si_per_mass(model, cross_section_2del_2del(model) / m_del(model));
  } else {
    set_rd_eta(model, NAN);
    set_rd_del(model, NAN);
    set_xi_cmb(model, NAN);
    set_xi_bbn(model, NAN);
    set_dneff_bbn(model, NAN);
    set_dneff_cmb(model, NAN);
    set_eta_si_per_mass(model, NAN);
    set_del_si_per_mass(model, NAN);
  }
}

} // namespace model

} // namespace darksun

#endif // DARKSUN_MODEL_BOLTZMANN_HPP
