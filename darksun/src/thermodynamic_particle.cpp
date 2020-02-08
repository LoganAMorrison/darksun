//
// Created by Logan Morrison on 2/7/20.
//

#include "darksun/thermodynamic_particle.hpp"
#include <boost/math/special_functions/bessel.hpp>
#include <cmath>

namespace darksun {


double ThermodynamicParticle::nbar(double x) const {
    if (x == 0.0)
        return (spin2 % 2 == 0) ? 0.121793828233573 : 0.0913453711751798;
    else {
        using boost::math::cyl_bessel_k;
        double eta = (spin2 % 2 == 0) ? 1.0 : -1.0;

        double sum = 0.0;
        for (int n = 1; n <= 5; n++) {
            sum += pow(eta, n + 1) * cyl_bessel_k(2, n * x) / n;
        }
        return x * x * sum / (2 * M_PI * M_PI);
    }
}

double ThermodynamicParticle::rhobar(double x) const {
    if (x == 0.0)
        return M_PI * M_PI / 30.0 * ((spin2 % 2 == 0) ? 1.0 : 7.0 / 8.0);
    else {
        using boost::math::cyl_bessel_k;
        double eta = (spin2 % 2 == 0) ? 1.0 : -1.0;

        double sum = 0.0;
        for (int n = 1; n <= 5; n++) {
            sum += pow(eta, n + 1) / (n * n) * (x * n * cyl_bessel_k(1, n * x) + 3.0 * cyl_bessel_k(2, n * x));
        }
        return x * x * sum / (2 * M_PI * M_PI);
    }
}

double ThermodynamicParticle::pbar(double x) const {
    if (x == 0.0)
        return M_PI * M_PI / 90.0 * ((spin2 % 2 == 0) ? 1 : 7.0 / 8.0);
    else {
        using boost::math::cyl_bessel_k;
        double eta = (spin2 % 2 == 0) ? 1.0 : -1.0;

        double sum = 0.0;
        for (int n = 1; n <= 5; n++) {
            sum += pow(eta, n + 1) / (n * n) * cyl_bessel_k(2, n * x);
        }
        return x * x * sum / (2 * M_PI * M_PI);
    }
}

double ThermodynamicParticle::sbar(double x) const {
    if (x == 0.0)
        return 2.0 * M_PI * M_PI / 45.0 * ((spin2 % 2 == 0) ? 1.0 : 7.0 / 8.0);
    else {
        using boost::math::cyl_bessel_k;
        double eta = (spin2 % 2 == 0) ? 1.0 : -1.0;

        double sum = 0.0;
        for (int n = 1; n <= 5; n++) {
            sum += pow(eta, n + 1) / n * cyl_bessel_k(3, n * x);
        }
        return x * x * x * sum / (2 * M_PI * M_PI);
    }
}

double ThermodynamicParticle::neq(double T) const {
    return g * nbar(mass / T) * T * T * T;
}

double ThermodynamicParticle::energy_density(double T) const {
    return g * rhobar(mass / T) * T * T * T * T;
}

double ThermodynamicParticle::pressure_density(double T) const {
    return g * pbar(mass / T) * T * T * T * T;
}

double ThermodynamicParticle::entropy_density(double T) const {
    return g * sbar(mass / T) * T * T * T;
}

double ThermodynamicParticle::dof_energy(double T) const {
    return 30.0 / (M_PI * M_PI) * g * rhobar(mass / T);
}

double ThermodynamicParticle::dof_entropy(double T) const {
    return g * 45.0 / (2.0 * M_PI * M_PI) * sbar(mass / T);
}

}
