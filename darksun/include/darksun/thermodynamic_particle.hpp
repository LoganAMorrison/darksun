/**
 * Created by Logan Morrison on 2019-06-09.
 *
 * Class for computing various thermodynamic quantities of a particle in
 * thermal equilibrium.
 *
 * Attributes:
 *      double mass: mass of the particle in GeV
 *      double g: internal dof of the particle
 *      double spin2: 2 times the particles spin (i.e., for spin=1/2, spin2 = 1)
 *
 * Methods:
 *      double neq(double T): equilibrium number density
 *      double energy_density(double T): equilibrium energy density
 *      double pressure_density(double T): equilibrium pressure density
 *      double entropy_density(double T): equilibrium entropy density
 *      double dof_energy(double T): equilibrium dof in energy
 *      double dof_entropy(double T): equilibrium dof in entropy
 *      
 */


#pragma once

namespace darksun {

class ThermodynamicParticle {
protected:
    /**
     * Number density scaled by the particle's temperature and internal d.o.f.
     * @param x mass divided by temperature: m / T
     * @return neq(T) / (g T^3)
     */
    double nbar(double x) const;

    /**
     * Energy density scaled by the particle's temperature and internal d.o.f.
     * @param x mass divided by temperature: m / T
     * @return energy_density(T) / (g T^4)
     */
    double rhobar(double x) const;

    /**
     * Pressure density scaled by the particle's temperature and internal d.o.f.
     * @param x mass divided by temperature: m / T
     * @return pressure_density(T) / (g T^4)
     */
    double pbar(double x) const;

    /**
     * Entropy density scaled by the particle's temperature and internal d.o.f.
     * @param x mass divided by temperature: m / T
     * @return entropy_density(T) / (g T^3)
     */
    double sbar(double x) const;

public:
    double mass;
    double g; // internal d.o.f.
    unsigned int spin2; // 2 * the particles spin

    ThermodynamicParticle(double mass, double g, unsigned int spin2) : mass(mass), g(g), spin2(spin2) {}

    ~ThermodynamicParticle() = default;

    /**
     * Compute the equilibrium number density of particle with temperature T.
     * @param T Temperature of particle
     * @return Number density
     */
    double neq(double T) const;

    /**
     * Compute the equilibrium energy density of particle with temperature T.
     * @param T Temperature of particle
     * @return Energy density
     */
    double energy_density(double T) const;

    /**
     * Compute the equilibrium pressure density of particle with temperature T.
     * @param T Temperature of particle
     * @return Pressure density
     */
    double pressure_density(double T) const;

    /**
     * Compute the equilibrium entropy density of particle with temperature T.
     * @param T Temperature of particle
     * @return Entropy density
     */
    double entropy_density(double T) const;

    /**
     * Compute the equilibrium number of d.o.f. stored in energy of a particle
     * with temperature T
     * @param T Temperature of particle
     * @return D.o.f. stored in energy
     */
    double dof_energy(double T) const;

    /**
     * Compute the equilibrium number of d.o.f. stored in entropy of a particle
     * with temperature T
     * @param T Temperature of particle
     * @return D.o.f. stored in entropy
     */
    double dof_entropy(double T) const;

};

}
