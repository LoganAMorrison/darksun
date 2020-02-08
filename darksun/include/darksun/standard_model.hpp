//
// Created by Logan Morrison on 1/26/20.
//

#pragma once

namespace darksun {

constexpr double _hsm0 = 3.9387999991430975;
constexpr double _hsm_inf = 106.83;

/**
 * Compute the effective number of degrees of freedom in energy of the SM.
 * @param T Temperature of the SM bath.
 * @return g_eff.
 */
double sm_dof_energy(double T);

/**
 * Compute the effective number of degrees of freedom in entropy of the SM.
 * @param T Temperature of the SM bath.
 * @return h_eff.
 */
double sm_dof_entropy(double T);

/**
 * Compute the value of h_eff / sqrt(g_eff) * (1 + (T / 3 h_eff) * dh_eff / dT)
 * @param Temperature of the SM bath.
 * @return Square-root of g_star.
 */
double sm_sqrt_gstar(double T);

/**
 * Compute the energy density of the SM bath.
 * @param T Temperature of the SM bath.
 * @return Energy density.
 */
double sm_energy_density(double T);

/**
 * Compute the entropy density of the SM bath.
 * @param T Temperature of the SM bath.
 * @return Entropy density.
 */
double sm_entropy_density(double T);

}
