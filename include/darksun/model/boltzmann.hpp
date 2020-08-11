#ifndef DARKSUN_MODEL_BOLTZMANN_HPP
#define DARKSUN_MODEL_BOLTZMANN_HPP

#include "darksun/constants.hpp"
#include "darksun/model/compute_xi.hpp"
#include "darksun/model/cross_sections.hpp"
#include "darksun/model/dneff.hpp"
#include "darksun/model/parameters.hpp"
#include "darksun/model/thermal_functions.hpp"
#include <fmt/core.h>
#include <gsl/gsl_errno.h>
#include <stiff/stiff.hpp>

namespace darksun {

//===========================================================================
//---- RHS of the Boltzmann -------------------------------------------------
//===========================================================================

void boltzmann(int *, double *t, double *y, double *dy,
               const DarkSunParameters &params) {

  const double x = exp(*t);
  const double meta = m_eta(params);
  const double tsm = meta / x;
  const double we = y[0]; // log(Yeta)

  const double xi = compute_xi_const_tsm(tsm, params);
  const double td = xi * tsm;
  const double s = StandardModel::entropy_density(tsm);
  const double we_eq = weq_eta(tsm, xi, params);

  const double sige = thermal_cross_section_4eta_2eta(meta / td, params);
  const double sigd = thermal_cross_section_2eta_2del(meta / td, params);

  const double com =
      sqrt(M_PI / 45) * M_PLANK * sqrt_gstar(tsm, xi, params) * tsm;
  const double pfe = -s * s * com;
  const double pfd = com / x;

  dy[0] = pfe * sige * exp(we) * (exp(2 * we) - exp(2 * we_eq));
  dy[1] = pfd * sigd * exp(2 * we);
}

//===========================================================================
//---- Jacobian RHS of the Boltzmann ----------------------------------------
//===========================================================================

void boltzmann_jac(int *, double *t, double *y, double *dfy, int *,
                   const DarkSunParameters &params) {

  const double x = exp(*t);
  const double meta = m_eta(params);
  const double tsm = meta / x;
  const double we = y[0]; // log(Yeta)

  const double xi = compute_xi_const_tsm(tsm, params);
  const double td = xi * tsm;
  const double s = StandardModel::entropy_density(tsm);
  const double we_eq = weq_eta(tsm, xi, params);

  const double sige = thermal_cross_section_4eta_2eta(meta / td, params);
  const double sigd = thermal_cross_section_2eta_2del(meta / td, params);

  const double com =
      sqrt(M_PI / 45) * M_PLANK * sqrt_gstar(tsm, xi, params) * tsm;
  const double pfe = com * (-s * s);
  const double pfd = com / x;

  // dfe / dWe
  dfy[0] = pfe * sige * exp(we) * (3.0 * exp(2 * we) - exp(2 * we_eq));
  // dfe / dYe
  dfy[1] = 0.0;
  // dfd / dWe
  dfy[2] = 2.0 * pfd * sigd * exp(2 * we);
  // dfd / dYd
  dfy[3] = 0.0;
}

//===========================================================================
//---- solution output ------------------------------------------------------
//===========================================================================

void solout(int *nr, double *logxold, double *logx, double *y, double *cont,
            int *lrc, int *, DarkSunParameters &params, int *,
            const stiff::RadauWeight &w) {

  // Determine if the eta' has frozen out
  const double x = exp(*logx);
  const double meta = m_eta(params);
  const double tsm = meta / x;
  const double we = y[0];
  const double xi = compute_xi_const_tsm(tsm, params);
  const double we_eq = weq_eta(tsm, xi, params);
  if (we - we_eq > 0.1 && params.xi_fo < 0.0) {
    params.xi_fo = xi;
    params.tsm_fo = tsm;
  }

  double dx = params.dlogx;
  double d = params.logx;
  int i = params.sol_idx;
  if (*nr == 1) {
    d = *logxold;
  }
  while ((*logxold <= d) && (*logx >= d)) {
    int idx1 = 1, idx2 = 2;
    params.ts[i] = d;
    params.ys[i][0] = contra(&idx1, &d, cont, lrc, w);
    params.ys[i][1] = contra(&idx2, &d, cont, lrc, w);

    d += dx;
    i += 1;
  }
  params.sol_idx = i;
  params.logx = d;
}

//===========================================================================
//---- Solve the Boltzmann --------------------------------------------------
//===========================================================================

void solve_boltzmann(double reltol, double abstol, DarkSunParameters &params) {
  using namespace stiff;

  // Initial conditions
  double meta = m_eta(params);

  // Integration interval
  double td = params.lam / 2.0;                // Start Td at confinement
  double xi = compute_xi_const_td(td, params); // Starting value of xi
  double tsm = td / xi;                        // Initial SM temperature
  double start = log(meta / tsm);
  double final = log(meta / T_CMB);
  params.dlogx = (final - start) / double(DarkSunParameters::SOL_LENGTH - 1);

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
  work[6] = 1e-2; // Set maximum step size to 0.01

  //==================================================================
  //---- Define lambdas which capture the model ----------------------
  //==================================================================
  auto boltz = [&params, &ipar](int *n, double *logx, double *y, double *dy) {
    boltzmann(n, logx, y, dy, params);
  };
  auto jac = [&params, &ipar](int *n, double *logx, double *y, double *dfy,
                              int *ldfy) {
    boltzmann_jac(n, logx, y, dfy, ldfy, params);
  };
  auto mas = [](int *, double *, int *) {};
  auto solo = [&params, &ipar](int *nr, double *xold, double *x, double *y,
                               double *cont, int *lrc, int *n, int *irtrn,
                               const RadauWeight &w) {
    solout(nr, xold, x, y, cont, lrc, n, params, irtrn, w);
  };

  //==================================================================
  //---- Solve the Boltzmann equations using RADAU -------------------
  //==================================================================
  radau(&nd, boltz, &start, y, &final, &h, &rtol, &atol, &itol, jac, &ijac,
        &mljac, &mujac, mas, &imas, &mlmas, &mumas, solo, &iout, work.data(),
        &lwork, iwork.data(), &liwork, &idid);

  //==================================================================
  //---- Record/Calculate outputs ------------------------------------
  //==================================================================
  // If 'idid' > 0, that means that RADAU was successful. Otherwise,
  // there was an error
  if (idid > 0) {
    // Save last state
    constexpr size_t last_idx = DarkSunParameters::SOL_LENGTH - 1;
    params.ts[last_idx] = final;
    params.ys[last_idx][0] = y[0];
    params.ys[last_idx][1] = y[1];

    params.rd_eta = m_eta(params) * exp(y[0]) * S_TODAY / RHO_CRIT;
    params.rd_del = m_del(params) * y[1] * S_TODAY / RHO_CRIT;
    params.xi_cmb = compute_xi_const_tsm(T_CMB, params);
    params.xi_bbn = compute_xi_const_tsm(T_BBN, params);
    params.dneff_cmb = compute_dneff_bbn(params);
    params.dneff_bbn = compute_dneff_bbn(params);
    params.eta_si_per_mass = cross_section_2eta_2eta(params) / m_eta(params);
    params.del_si_per_mass = cross_section_2del_2del(params) / m_del(params);
  } else {
    params.xi_fo = NAN;
    params.tsm_fo = NAN;
    params.rd_eta = NAN;
    params.rd_del = NAN;
    params.xi_cmb = NAN;
    params.xi_bbn = NAN;
    params.dneff_cmb = NAN;
    params.dneff_bbn = NAN;
    params.eta_si_per_mass = NAN;
    params.del_si_per_mass = NAN;
  }
}

} // namespace darksun

#endif // DARKSUN_MODEL_BOLTZMANN_HPP
