//
// Created by Logan Morrison on 2/5/20.
//
//#include "darksun/c"
#include "darksun/standard_model.hpp"
#include "darksun/thermodynamic_particle.hpp"
#include "darksun/odeint.hpp"
//#include "darksun/odeint_radau.hpp"
//#include "darksun/stepper_dopr853.hpp"
#include "darksun/stepper_sie.hpp"
#include <boost/math/differentiation/finite_difference.hpp>
#include <Eigen/Dense>
#include <chrono>
#include <iostream>

using namespace std::chrono;
using namespace darksun;

constexpr double PLANK_MASS = 1.220910e19;
constexpr double CRITICAL_ENERGY_DENSITY = 1.05375e-5;
constexpr double SM_ENTROPY_DENSITY_TODAY = 2891.2;

class BoltzmannTest {
private:
public:
    ThermodynamicParticle chi;
    double mdm;
    double sigmav;

    BoltzmannTest(double mdm, double sigmav) : chi(ThermodynamicParticle{mdm, 1.0, 1}), mdm(mdm), sigmav(sigmav) {}

    ~BoltzmannTest() = default;

    double operator()(double logx, Eigen::VectorXd &w, Eigen::VectorXd &dwdlogx) {
        double x = exp(logx);
        double T = mdm / x;
        double pf = -sqrt(M_PI / 45) * PLANK_MASS * mdm * sm_sqrt_gstar(T) / x;
        double weq = log(chi.neq(T) / sm_entropy_density(T));

        dwdlogx(0) = pf * sigmav * (exp(w(0)) - exp(2.0 * weq - w(0)));
        return dwdlogx(0);
    }

    void jacobian(double logx, Eigen::VectorXd &w, Eigen::VectorXd &dfdx, Eigen::MatrixXd &dfdy) {
        using namespace boost::math::differentiation;
        double x = exp(logx);
        double T = mdm / x;
        double pf = -sqrt(M_PI / 45) * PLANK_MASS * mdm * sm_sqrt_gstar(T) / x;
        double weq = log(chi.neq(T) / sm_entropy_density(T));

        dfdx(0) = finite_difference_derivative(
                [this, &w, dfdy](double logx) {
                    Eigen::VectorXd dwdlogx(1);
                    return operator()(logx, w, dwdlogx);
                }, logx);
        dfdy(0, 0) = pf * sigmav * (exp(w(0)) + exp(2.0 * weq - w(0)));
    }
};

int main() {
    /*
     * mx1, sigmav1 = 10.313897683787216e3, 1.966877938634266e-15
     * mx2, sigmav2 = 104.74522360006331e3, 1.7597967261428258e-15
     * mx3, sigmav3 = 1063.764854316313e3, 1.837766552668581e-15
     * mx4, sigmav4 = 10000.0e3, 1.8795945459427076e-15
     */
    //BoltzmannTest model{10.313897683787216, 1.966877938634266e-9};
    BoltzmannTest model{10000.0, 1.8795945459427076e-9};

    double x_start = 1.0;
    Eigen::VectorXd wstart(1);
    wstart(0) = log(model.chi.neq(model.mdm / x_start) / sm_entropy_density(model.mdm / x_start));

    std::cout << "w(logx=" << log(x_start) << ") = " << wstart(0) << std::endl;


    Output output(100);
    Odeint<StepperSie<BoltzmannTest>> solver(wstart, log(x_start), log(x_start) + 7.0, 1e-5, 1e-3, 1e-3, 0.0, output,
                                             model);

    auto start = high_resolution_clock::now();

    solver.integrate();

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<nanoseconds>(stop - start);
    std::cout << "duration: = " << duration.count() << std::endl;

    for (int i = 0; i < output.count; i++) {
        std::cout << "logx = " << output.xsave[i] << ", w = " << output.ysave(0, i) << std::endl;
    }

    double yinf = exp(output.ysave(0, output.count - 1));

    double rd = yinf * model.mdm * SM_ENTROPY_DENSITY_TODAY / CRITICAL_ENERGY_DENSITY;

    std::cout << std::endl;
    std::cout << "winf = " << output.ysave(0, output.count - 1) << std::endl;
    std::cout << "rd = " << rd << std::endl;

    return 0;
}

