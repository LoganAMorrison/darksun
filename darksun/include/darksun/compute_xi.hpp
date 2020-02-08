//
// Created by Logan Morrison on 1/27/20.
//

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
protected:
private:
public:

    double xi_inf = 1.0; // ratio of dark to SM temperature at above EW scale
    double hd_inf = 1.0; // value of entropic dof in dark sector above EW scale
    double sum_g = 1.0; // sum of internal dof in the dark sector
    double gl = 1.0; // internal dof of the lightest particle in dark sector
    double ml = 0.0; // mass of the lightest dof in dark sector

    ThermallyDecoupledModel() = default;

    ~ThermallyDecoupledModel() = default;

    virtual double dark_dof_entropy(double T) const = 0;

    double xi_upper_bound_const_td(double Td) const;

    double xi_lower_bound_const_td(double Td) const;

    double xi_upper_bound_const_tsm(double Tsm) const;

    double xi_lower_bound_const_tsm(double Tsm) const;

    double compute_xi_const_td(double Td) const;

    double compute_xi_const_tsm(double Tsm) const;
};

}
