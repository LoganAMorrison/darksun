#ifndef DARK_SUN_PHASE_SPACE_RAMBO_HPP
#define DARK_SUN_PHASE_SPACE_RAMBO_HPP

#include "darksun/phase_space/base.hpp"
#include "darksun/phase_space/four_momentum.hpp"
#include <functional>
#include <vector>

namespace darksun {

class Rambo : public PhaseSpaceGenerator {
private:
  double compute_scale_factor(std::vector<FourMomentum> &);

  void initialize_four_momenta(std::vector<FourMomentum> &);

  void boost_four_momenta(std::vector<FourMomentum> &);

  double correct_masses(std::vector<FourMomentum> &);

  PhaseSpaceEvent internal_generate_event();

  void internal_generate_events(std::size_t);

public:
  Rambo(std::vector<double> &isp_masses, std::vector<double> &fsp_masses,
        double cme);

  Rambo(std::vector<double> &isp_masses, std::vector<double> &fsp_masses,
        double cme,
        std::function<double(const std::vector<FourMomentum> &)> t_mat_squared);

  PhaseSpaceEvent generate_event();

  std::vector<PhaseSpaceEvent> generate_events(std::size_t) override;
};

Rambo::Rambo(std::vector<double> &isp_masses, std::vector<double> &fsp_masses,
             double cme)
    : PhaseSpaceGenerator(isp_masses, fsp_masses, cme) {}

Rambo::Rambo(
    std::vector<double> &isp_masses, std::vector<double> &fsp_masses,
    double cme,
    std::function<double(const std::vector<FourMomentum> &)> t_mat_squared)
    : PhaseSpaceGenerator(isp_masses, fsp_masses, cme,
                          std::move(t_mat_squared)) {}

/* Function for finding the scaling parameter to turn mass-less four-vectors
 * into four-vectors with the correct masses.
 * @param momenta 4-momenta of final-state particles
 */
double Rambo::compute_scale_factor(std::vector<FourMomentum> &momenta) {
  static thread_local const int MAX_ITER = 50;
  static thread_local const double TOL = 1e-4;

  double mass_sum = std::accumulate(fsp_masses.begin(), fsp_masses.end(), 0.0);

  double xi = sqrt(1.0 - (mass_sum / cme) * (mass_sum / cme));

  int iter_count = 0;
  bool converged = false;
  do { // Perform newton iterations to solve for xi
    double f = -cme;
    double df = 0.0;

    for (size_t i = 0; i < fsp_masses.size(); i++) {
      // Compute residual and derivative of residual
      double m2 = fsp_masses[i] * fsp_masses[i];
      double xi2 = xi * xi;
      double e2 = momenta[i].e * momenta[i].e;
      double del_f = sqrt(m2 + xi2 * e2);
      f += del_f;
      df += xi * e2 / del_f;
    }

    // Newton correction
    double delta_xi = -(f / df);
    xi += delta_xi;

    iter_count++;
    if (fabs(delta_xi) < TOL || iter_count >= MAX_ITER) {
      converged = true;
    }
  } while (!converged);
  return xi;
}

/**
 * Initialize the four-momenta with isotropic, random four-momenta with
 * energies, q₀, distributed according to q₀ * exp(-q₀).
 * @param momenta 4-momenta of final-state particles
 */
void Rambo::initialize_four_momenta(std::vector<FourMomentum> &momenta) {

  std::uniform_real_distribution<double> distribution(0.0, 1.0);

  for (size_t i = 0; i < fsp_masses.size(); i++) {
    double rho1 = phase_space_uniform_rand();
    double rho2 = phase_space_uniform_rand();
    double rho3 = phase_space_uniform_rand();
    double rho4 = phase_space_uniform_rand();

    double c = 2.0 * rho1 - 1.0;
    double phi = 2.0 * M_PI * rho2;

    momenta[i].e = -log(rho3 * rho4);
    momenta[i].px = momenta[i].e * sqrt(1.0 - c * c) * cos(phi);
    momenta[i].py = momenta[i].e * sqrt(1.0 - c * c) * sin(phi);
    momenta[i].pz = momenta[i].e * c;
  }
}

/**
 * Boost the four-momenta into the center-of-mass frame and compute the
 * initial weight of the event.
 * @param momenta 4-momenta of final-state particles
 */
void Rambo::boost_four_momenta(std::vector<FourMomentum> &momenta) {
  // Total momentum and its mass
  FourMomentum Q =
      std::accumulate(momenta.begin(), momenta.end(), FourMomentum{});
  double massQ = Q.mass();

  // Boost three-vector
  double bx = -Q.px / massQ;
  double by = -Q.py / massQ;
  double bz = -Q.pz / massQ;
  // Boost factors
  double x = cme / massQ;
  double gamma = Q.e / massQ;
  double a = 1.0 / (1.0 + gamma);

  for (size_t i = 0; i < fsp_masses.size(); i++) {
    double qe = momenta[i].e;
    double qx = momenta[i].px;
    double qy = momenta[i].py;
    double qz = momenta[i].pz;

    double b_dot_q = bx * qx + by * qy + bz * qz;

    momenta[i].e = x * (gamma * qe + b_dot_q);
    momenta[i].px = x * (qx + bx * qe + a * b_dot_q * bx);
    momenta[i].py = x * (qy + by * qe + a * b_dot_q * by);
    momenta[i].pz = x * (qz + bz * qe + a * b_dot_q * bz);
  }
}

/**
 * Correct the masses of the four-momenta and correct the weight of the
 * event.
 * @param momenta 4-momenta of final-state particles
 * @return new event weight factor
 */
double Rambo::correct_masses(std::vector<FourMomentum> &momenta) {
  double xi = compute_scale_factor(momenta);

  double term1 = 0.0;
  double term2 = 0.0;
  double term3 = 1.0;

  for (size_t i = 0; i < fsp_masses.size(); i++) {
    double m = fsp_masses[i];
    double eng = momenta[i].e;
    momenta[i].e = sqrt(m * m + (xi * eng) * (xi * eng));
    momenta[i].px *= xi;
    momenta[i].py *= xi;
    momenta[i].pz *= xi;

    double mod =
        sqrt(momenta[i].px * momenta[i].px + momenta[i].py * momenta[i].py +
             momenta[i].pz * momenta[i].pz);
    eng = momenta[i].e;

    term1 += mod / cme;
    term2 += mod * mod / eng;
    term3 *= mod / eng;
  }

  term1 = pow(term1, 2.0 * fsp_masses.size() - 3.0);
  term2 = 1.0 / term2;

  // re-weight
  return term1 * term2 * term3 * cme;
}

/**
 * generate single phase space event
 * @return event
 */
PhaseSpaceEvent Rambo::internal_generate_event() {
  std::vector<FourMomentum> momenta(fsp_masses.size(), FourMomentum{});
  double weight;

  initialize_four_momenta(momenta);
  boost_four_momenta(momenta);
  weight = correct_masses(momenta) * mat_squared(momenta) * m_base_weight;

  return PhaseSpaceEvent{momenta, weight};
}

/**
 * Generate many phase space events.
 * @param num_points
 * @return vector of events
 */
void Rambo::internal_generate_events(size_t num_points) {
  std::vector<PhaseSpaceEvent> local_events;
  local_events.reserve(num_points);

  for (size_t n = 0; n < num_points; n++)
    local_events.emplace_back(internal_generate_event());

  {
    // Add events to class level event array. Avoid data races by lock guard.
    std::lock_guard<std::mutex> lck(m_mtx);
    for (auto &event : local_events)
      m_events.emplace_back(event);
  }
}

/**
 * Generate a set of Rambo event.
 * @return RamboEvent.
 */
PhaseSpaceEvent Rambo::generate_event() {
  auto num_fsp_d = (double)fsp_masses.size();
  m_base_weight = pow(M_PI / 2.0, num_fsp_d - 1.0) *
                  pow(cme, 2.0 * num_fsp_d - 4.0) / tgamma(num_fsp_d) /
                  tgamma(num_fsp_d - 1.0) *
                  pow(2.0 * M_PI, 4.0 - 3.0 * num_fsp_d);
  return internal_generate_event();
}

/**
 * Generate a set of Rambo events.
 * @param num_events number of events to generate.
 * @return nothing; events stored in 'events'.
 */
std::vector<PhaseSpaceEvent> Rambo::generate_events(size_t num_events) {
  auto num_fsp_d = (double)fsp_masses.size();
  m_base_weight = pow(M_PI / 2.0, num_fsp_d - 1.0) *
                  pow(cme, 2.0 * num_fsp_d - 4.0) / tgamma(num_fsp_d) /
                  tgamma(num_fsp_d - 1.0) *
                  pow(2.0 * M_PI, 4.0 - 3.0 * num_fsp_d);

  m_events.clear();
  size_t num_threads = std::thread::hardware_concurrency();
  std::vector<std::thread> threads;

  for (size_t n = 0; n < num_threads - 1; n++) {
    threads.emplace_back(
        [this](size_t num_points) {
          this->internal_generate_events(num_points);
        },
        num_events / num_threads);
  }
  // add the remainder of points into last thread
  threads.emplace_back(
      [this](size_t num_points) { this->internal_generate_events(num_points); },
      num_events / num_threads + (num_events % num_threads));
  for (auto &thread : threads) {
    thread.join();
  }

  return m_events;
}

} // namespace darksun

#endif // DARK_SUN_PHASE_SPACE_RAMBO_HPP
