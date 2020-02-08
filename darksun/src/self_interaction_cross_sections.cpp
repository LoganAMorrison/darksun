//
// Created by Logan Morrison on 2/7/20.
//

#include "darksun/models.hpp"
#include <cmath>
#include <boost/math/special_functions/pow.hpp>

namespace darksun {

double DarkSUNModel::cross_section_2eta_2eta() const {
    using namespace boost::math;

    //compute dark temperature today in order to compute eta velocity
    double Tsm_today = 2.7255 * 8.6173303e-14;
    double xi_today = compute_xi(Tsm_today);
    double Td_today = xi_today * Tsm_today;

    double me = eta.mass;
    double vrel = sqrt(3.0 * Td_today / me);
    double s = 4 * me * me / (1 - vrel * vrel);

    double fe = lam * sqrt(N) / (4.0 * M_PI);
    return (8.0 * L1 * L1 * (376.0 * pow<8>(me) - 576.0 * pow<6>(me) * s + 396.0 * pow<4>(me) * s * s
            - 136 * me * me * pow<3>(s) + 21 * pow<4>(s))) / (15.0 * pow<4>(fe * lam));
}

double DarkSUNModel::cross_section_2delta_2delta() const {
    /* diagrams with sigma exchange
     * ----------------------------
     * p1 -->--.-->-- p3       p1 -->--.-->-- p4
     *         .                       .
     *         .          +            .
     *         .                       .
     * p2 -->--.-->-- p4       p2 -->--.-->-- p3
     */
    using namespace boost::math;
    double g = 1.0; // TODO: not sure what to use for this
    double msigma = lam / sqrt(N);
    double mdelta = delta.mass;

    return 1.5 * pow<4>(g) * N * mdelta * mdelta / (32.0 * M_PI) / pow<4>(msigma);
}


}

