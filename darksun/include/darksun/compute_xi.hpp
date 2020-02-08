/* Created by Logan Morrison on 1/27/20.
 *
 * The "c" abstract base class defines functions to compute the
 * temperature ratio between the standard model bath, and a secluded (decoupled)
 * sector. The temperature ratio, xi, is defined as the ratio of decoupled to
 * standard model temperature: xi = T_d / T_sm.
 *
 * Classes deriving from ThermallyDecoupledModel must define the following:
 *      double xi_inf: ratio of dark to SM temperature at above EW scale
 *      double hd_inf: value of entropic dof in dark sector above EW scale
 *      double sum_g: sum of internal dof in the dark sector
 *      double gl: internal dof of the lightest particle in dark sector
 *      double ml: mass of the lightest dof in dark sector
 *
 *      double dark_dof_entropy(double Td): function to compute the entropic
 *                                          dof in the dark sector
 *
 * The temperature ratio, xi = Td / Tsm, can be computed in two ways. The
 * first, if T_d is known, then we can determine xi using:
 *      compute_xi_const_td(Td)
 * second, if T_sm is known, then:
 *      compute_xi_const_tsm(Tsm)
 */


#pragma once

#include <cmath>
#include <functional>
#include <boost/math/special_functions/lambert_w.hpp>
#include <boost/math/tools/roots.hpp>
#include "darksun/standard_model.hpp"

namespace darksun {

/**
 * Abstract base class for computing the ratio of dark to SM temperature for a
 * model in which the dark sector is thermally decoupled from the SM.
 */
class ThermallyDecoupledModel {
private:

    /**
     * Compute the upper bound on xi assuming Td is known.
     * @param Td Temperature of the dark bath
     * @return Upper bound on xi
     */
    double xi_upper_bound_const_td(double Td) const;

    /**
     * Compute the lower bound on xi assuming Td is known.
     * @param Td Temperature of the dark bath
     * @return Lower bound on xi
     */
    double xi_lower_bound_const_td(double Td) const;

    /**
     * Compute the upper bound on xi assuming Tsm is known.
     * @param Tsm Temperature of the standard model bath
     * @return Upper bound on xi
     */
    double xi_upper_bound_const_tsm(double Tsm) const;

    /**
     * Compute the lower bound on xi assuming Tsm is known.
     * @param Tsm Temperature of the standard model bath
     * @return Lower bound on xi
     */
    double xi_lower_bound_const_tsm(double Tsm) const;
public:

    double xi_inf = 1.0; // ratio of dark to SM temperature at above EW scale
    double hd_inf = 1.0; // value of entropic dof in dark sector above EW scale
    double sum_g = 1.0; // sum of internal dof in the dark sector
    double gl = 1.0; // internal dof of the lightest particle in dark sector
    double ml = 0.0; // mass of the lightest dof in dark sector

    ThermallyDecoupledModel() = default;

    ~ThermallyDecoupledModel() = default;

    /**
    * Compute the entropic dof of the dark sector
    * @param Td Temperature of the dark bath
    * @return h_eff entropic dof of the dark sector
    */
    virtual double dark_dof_entropy(double Td) const = 0;

    /**
     * Compute xi = Td/Tsm assuming Td is known
     * @param Td Temperature of the dark bath
     * @return Value of xi = Td / Tsm
     */
    double compute_xi_const_td(double Td) const;

    /**
     * Compute xi assuming Tsm is known
     * @param Tsm Temperature of the standard model bath
     * @return Value of xi = Td / Tsm
     */
    double compute_xi_const_tsm(double Tsm) const;
};

}
