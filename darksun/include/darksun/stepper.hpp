//
// Created by Logan Morrison on 2/5/20.
//

#pragma once

#include <Eigen/Dense>

namespace darksun {
struct StepperBase {
    double &x;
    double xold{};
    Eigen::VectorXd &y, &dydx;
    double atol, rtol;
    bool dense;
    double hdid{};
    double hnext{};
    double EPS{};
    int n, neqn;
    Eigen::VectorXd yout, yerr;

    StepperBase(Eigen::VectorXd &yy, Eigen::VectorXd &dydxx, double &xx, const double atoll,
                const double rtoll, bool dens)
            : x(xx), y(yy), dydx(dydxx), atol(atoll),
              rtol(rtoll), dense(dens), n(y.size()), neqn(n), yout(n), yerr(n) {}
};

}
