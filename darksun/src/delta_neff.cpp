//
// Created by Logan Morrison on 2/7/20.
//

#include "darksun/models.hpp"
#include "darksun/standard_model.hpp"
#include <boost/math/special_functions/pow.hpp>

namespace darksun {

double DarkSUNModel::delta_n_eff_cmb() const {
    using namespace boost::math;
    return 4.0 / 7.0 * pow(11.0 / 4.0, 4.0 / 3.0) *
            dark_dof_energy(T_CMB * xi_cmb) * pow<4>(xi_cmb);
}

double DarkSUNModel::delta_n_eff_bbn() const {
    using namespace boost::math;
    return 4.0 / 7.0 * dark_dof_energy(T_BBN * xi_bbn) * pow<4>(xi_bbn);
}

}
