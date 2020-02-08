//
// Created by Logan Morrison on 2019-05-24.
//

#include "darksun/thermodynamic_particle.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <cmath>
#include <boost/math/quadrature/gauss_kronrod.hpp>

using namespace darksun;

constexpr double log_x_min = -2.0;
constexpr double log_x_max = 2.0;
constexpr int num_xs = 100;

double neq_exact(double T, ThermodynamicParticle &p);

double energy_density_exact(double T, ThermodynamicParticle &p);

double pressure_density_exact(double T, ThermodynamicParticle &p);

double entropy_density_exact(double T, ThermodynamicParticle &p);

/**
 * Class for testing the Fermion class. We used Mathematica to compute the
 * data. All the data has the temperature and d.o.f. terms factored out. For
 * example:
 *      nbar = number_density(T) / (g T^3)
 */
class FermionTest : public ::testing::Test {
protected:
    void SetUp() override {
        xs.reserve(num_xs);
        double step_log_x = (log_x_max - log_x_min) / (num_xs - 1);
        for (int n = 0; n < num_xs; n++) {
            double log_x = log_x_min + n * step_log_x;
            xs.push_back(pow(10.0, log_x));
        }
    }

    // void TearDown() override {}
    std::vector<double> xs;
    ThermodynamicParticle fermion{1.0, 1.0, 1};
};

/**
 * Class for testing the Boson class. We used Mathematica to compute the
 * data. All the data has the temperature and d.o.f. terms factored out. For
 * example:
 *      nbar = number_density(T) / (g T^3)
 */
class BosonTest : public ::testing::Test {
protected:
    void SetUp() override {
        xs.reserve(num_xs);
        double step_log_x = (log_x_max - log_x_min) / (num_xs - 1);
        for (int n = 0; n < num_xs; n++) {
            double log_x = log_x_min + n * step_log_x;
            xs.push_back(pow(10.0, log_x));
        }
    }

    // void TearDown() override {}
    std::vector<double> xs;
    ThermodynamicParticle boson{1.0, 1.0, 0};
};

/**
 * Test that the Fermion class correctly compute the equilibrium number density.
 */
TEST_F(FermionTest, TestNumberDensity) {
    for (double x : xs) {
        double T = fermion.mass / x;
        double approx = fermion.neq(T);
        double exact = neq_exact(T, fermion);
        double frac_diff = fabs(approx - exact) / exact;
        ASSERT_LE(frac_diff, 5e-2);
    }
}

/**
 * Test that the Fermion class correctly compute the equilibrium energy density.
 */
TEST_F(FermionTest, TestEnergyDensity) {
    for (double x : xs) {
        double T = fermion.mass / x;
        double approx = fermion.energy_density(T);
        double exact = energy_density_exact(T, fermion);
        double frac_diff = fabs(approx - exact) / exact;
        ASSERT_LE(frac_diff, 5e-2);
    }
}

/**
 * Test that the Fermion class correctly compute the equilibrium pressure density.
 */
TEST_F(FermionTest, TestPressureDensity) {
    for (double x : xs) {
        double T = fermion.mass / x;
        double approx = fermion.pressure_density(T);
        double exact = pressure_density_exact(T, fermion);
        double frac_diff = fabs(approx - exact) / exact;
        ASSERT_LE(frac_diff, 5e-2);
    }
}

/**
 * Test that the Fermion class correctly compute the equilibrium entropy density.
 */
TEST_F(FermionTest, TestEntropyDensity) {
    for (double x : xs) {
        double T = fermion.mass / x;
        double approx = fermion.entropy_density(T);
        double exact = entropy_density_exact(T, fermion);
        double frac_diff = fabs(approx - exact) / exact;
        ASSERT_LE(frac_diff, 5e-2);
    }
}

/**
 * Test that the Fermion class correctly compute the d.o.f. stored in energy.
 */
TEST_F(FermionTest, TestDOFEnergy) {
    for (double x : xs) {
        double T = fermion.mass / x;
        double approx = fermion.dof_energy(T);
        double exact = 30.0 / (M_PI * M_PI) * energy_density_exact(T, fermion) / pow(T, 4);
        double frac_diff = fabs(approx - exact) / exact;
        ASSERT_LE(frac_diff, 5e-2);
    }
}

/**
 * Test that the Fermion class correctly compute the d.o.f. stored in entropy.
 */
TEST_F(FermionTest, TestDOFEntropy) {
    for (double x : xs) {
        double T = fermion.mass / x;
        double approx = fermion.dof_energy(T);
        double exact = 45.0 / (2.0 * M_PI * M_PI) * entropy_density_exact(T, fermion) / pow(T, 3);
        double frac_diff = fabs(approx - exact) / exact;
        ASSERT_LE(frac_diff, 20e-2);
    }
}

/**
 * Test that the Boson class correctly compute the equilibrium number density.
 */
TEST_F(BosonTest, TestNumberDensity) {
    for (double x : xs) {
        double T = boson.mass / x;
        double approx = boson.neq(T);
        double exact = neq_exact(T, boson);
        double frac_diff = fabs(approx - exact) / exact;
        ASSERT_LE(frac_diff, 5e-2);
    }
}

/**
 * Test that the Boson class correctly compute the equilibrium energy density.
 */
TEST_F(BosonTest, TestEnergyDensity) {
    for (double x : xs) {
        double T = boson.mass / x;
        double approx = boson.energy_density(T);
        double exact = energy_density_exact(T, boson);
        double frac_diff = fabs(approx - exact) / exact;
        ASSERT_LE(frac_diff, 5e-2);
    }
}

/**
 * Test that the Boson class correctly compute the equilibrium pressure density.
 */
TEST_F(BosonTest, TestPressureDensity) {
    for (double x : xs) {
        double T = boson.mass / x;
        double approx = boson.pressure_density(T);
        double exact = pressure_density_exact(T, boson);
        double frac_diff = fabs(approx - exact) / exact;
        ASSERT_LE(frac_diff, 5e-2);
    }
}

/**
 * Test that the Fermion class correctly compute the equilibrium entropy density.
 */
TEST_F(BosonTest, TestEntropyDensity) {
    for (double x : xs) {
        double T = boson.mass / x;
        double approx = boson.entropy_density(T);
        double exact = entropy_density_exact(T, boson);
        double frac_diff = fabs(approx - exact) / exact;
        ASSERT_LE(frac_diff, 5e-2);
    }
}

/**
 * Test that the Boson class correctly compute the d.o.f. stored in energy.
 */
TEST_F(BosonTest, TestDOFEnergy) {
    for (double x : xs) {
        double T = boson.mass / x;
        double approx = boson.dof_energy(T);
        double exact = 30.0 / (M_PI * M_PI) * energy_density_exact(T, boson) / pow(T, 4);
        double frac_diff = fabs(approx - exact) / exact;
        ASSERT_LE(frac_diff, 5e-2);
    }
}

/**
 * Test that the Boson class correctly compute the d.o.f. stored in entropy.
 */
TEST_F(BosonTest, TestDOFEntropy) {
    for (double x : xs) {
        double T = boson.mass / x;
        double approx = boson.dof_energy(T);
        double exact = 45.0 / (2.0 * M_PI * M_PI) * entropy_density_exact(T, boson) / pow(T, 3);
        double frac_diff = fabs(approx - exact) / exact;
        ASSERT_LE(frac_diff, 10e-2);
    }
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}


double neq_exact(double T, ThermodynamicParticle &p) {
    using boost::math::quadrature::gauss_kronrod;

    double eta = (p.spin2 % 2 == 0) ? 1.0 : -1.0;

    auto integrand = [&p, eta, T](double eng) {
        return eng * sqrt(eng * eng - p.mass * p.mass) / (exp(eng / T) - eta);
    };

    double integral = gauss_kronrod<double, 15>::integrate(
            integrand, p.mass, std::numeric_limits<double>::infinity(), 5, 1e-9);

    return p.g * integral / (2.0 * M_PI * M_PI);
}

double energy_density_exact(double T, ThermodynamicParticle &p) {
    using boost::math::quadrature::gauss_kronrod;

    double eta = (p.spin2 % 2 == 0) ? 1.0 : -1.0;

    auto integrand = [&p, eta, T](double eng) {
        return eng * eng * sqrt(eng * eng - p.mass * p.mass) / (exp(eng / T) - eta);
    };

    double integral = gauss_kronrod<double, 15>::integrate(
            integrand, p.mass, std::numeric_limits<double>::infinity(), 5, 1e-9);

    return p.g * integral / (2.0 * M_PI * M_PI);
}

double pressure_density_exact(double T, ThermodynamicParticle &p) {
    using boost::math::quadrature::gauss_kronrod;

    double eta = (p.spin2 % 2 == 0) ? 1.0 : -1.0;

    auto integrand = [&p, eta, T](double eng) {
        return pow(eng * eng - p.mass * p.mass, 1.5) / (exp(eng / T) - eta);
    };

    double integral = gauss_kronrod<double, 15>::integrate(
            integrand, p.mass, std::numeric_limits<double>::infinity(), 5, 1e-9);

    return p.g * integral / (6.0 * M_PI * M_PI);
}

double entropy_density_exact(double T, ThermodynamicParticle &p) {
    return (energy_density_exact(T, p) + pressure_density_exact(T, p)) / T;
}