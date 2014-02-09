#ifndef LENNARD_JONES_SYSTEM_H_ZDNDUQM8

#define LENNARD_JONES_SYSTEM_H_ZDNDUQM8



// C++ Standard Library
#include <vector>

// Boost Array
#include <boost/array.hpp>

// Boost Program Options
#include <boost/program_options.hpp>

// Boost Random
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_on_sphere.hpp>
#include <boost/random/uniform_smallint.hpp>


#include "simulation_system.h"
#include "q_matrix_tools.h"
#include "rhab/pbc.h"


using namespace rhab;


class LennardJonesSystem : public SimulationSystem {
private:
  /**
   * Settings
   */
  // the sampler chosen
  size_t sampler;

  // size L of one dimension of the box
  double size;

  // number of particles
  size_t num_particles;

  // value sigma of Lennard-Jones potential
  double sigma;

  // value epsilon of Lennard-Jones potential
  double epsilon;

  // cufoff-radius
  double cutoff_radius;

  // initial delta r
  double delta_r;

  /**
   * Precalculated constants
   */
  // 4*epsilon
  double four_epsilon;

  // sigma**6
  double sigma6;

  // cutoff squared
  double cutoff_radius_sq;

  // dimensions of the box
  boost::array<double, 3> box_dimensions;

  /**
   * State
   */
  std::vector< boost::array<double, 3> > particles;

  /**
   * Misc
   */
  boost::uniform_01<> dist01;
  boost::uniform_smallint<size_t> ran_particle;
  typedef boost::uniform_on_sphere<double, boost::numeric::ublas::vector<double> > u_sphere_dist_t;
  u_sphere_dist_t u_sphere_dist;

  // exact dos
  vector_double_t dos_exact_norm;

  // size dependent constants:
  double e_min;
  double e_max;

  void set_bins(size_t n_bins);
  void set_particles(size_t num_particles);

  double calculate_dist_sq(const boost::array<double, 3> &a, const boost::array<double, 3> &b) {
    double r = 0;
    for (size_t i = 0; i < a.size(); ++i) {
      double tmp = a[i] - b[i];
      r += tmp*tmp;
    }
    return r;
  }

  double calculate_energy() {
    double energy(0.0);
    for (size_t i = 0; i < particles.size()-1; ++i) {
      for (size_t j = i+1; j < particles.size(); ++j) {
        double rsq = calculate_dist_sq(particles[i], particles[j]);
        energy += calculate_lj_interaction_rsq(rsq);
      }
    }
    return energy;
  }

  double calculate_lj_interaction(const double &r) {
    return calculate_lj_interaction_rsq(r*r);
  }

  double calculate_lj_interaction_rsq(const double &r2) {
    double r6 = r2*r2*r2;
    return four_epsilon * sigma6/r6 * (sigma6/r6 - 1);
  }

  template<class Sampler>
  void mc_loop() {
    for (size_t run = 0; run < runs; run++) {
      // variables for error / statistics calculation
      size_t error_check_freq = error_check_f;
      size_t index = 0;
      // reset Q matrix
      Q *= 0;
      double energy = calculate_energy();
      //int i, j;

      Sampler sampler(rng.rng, Q, settings);
      for (size_t step = 1; step <= steps; step++) {
        // select a particle at random
        size_t i = ran_particle(rng.rng);

        // offset
        u_sphere_dist_t::result_type offset = delta_r * u_sphere_dist(rng.rng);

        // calculate change in energy
        int dE = particles[i];

        int Ei = (energy-e_min) / 4;
        int Ej = (energy-e_min+dE) / 4;

        bool out_of_bounds = false;
        // update Q matrix
        if (Ej < 0 || Ej >= (int)n_bins) {
          Q(Ei, Ei)++;
          out_of_bounds = true;
          sampler.rejected(Ei);
        } else {
          Q(Ei, Ej)++;
        }

        if (!out_of_bounds && sampler(dE, Ei, Ej)) {
          //lattice(i, j) *= -1;
          energy += dE;
        }
        if (step % error_check_freq == 0) {
          error_acc[index].step = step;

          double err_lq, err_gth, err_power;
          bool   suc_lq, suc_gth, suc_power;
          double q_error;
          boost::tie(err_lq, err_gth, err_power, suc_lq, suc_gth, suc_power, q_error) = rhab::calculate_error_q(dos_exact_norm, matrix_double_t(), Q, Qd, error_matrices, index);
          if (suc_lq) {
            error_acc[index].err1(err_lq);
          }
          if (suc_gth) {
            error_acc[index].err2(err_gth);
          }
          if (suc_power) {
            error_acc[index].err3(err_power);
          }
          if (sampler.has_own_statistics()) {
            double err = sampler.calculate_error(dos_exact_norm, error_matrices, index);
            error_acc[index].err4(err);
          }

          index++;
          if (step % (10*error_check_freq) == 0) {
            error_check_freq *= 10;
          }
        }
        sampler.check(step, run);
      }
      //std::cout << Q << std::endl;
      //std::cout << Qd << std::endl;
      //std::cout << calculate_dos_power(Qd) << std::endl;
      //std::cout << dos_exact_norm << std::endl;
    }
  }

  //void read_exact_dos();

  void setup_output();

public:
  LennardJonesSystem(settings_t & s);

  virtual boost::program_options::options_description get_program_options();
  virtual void parse_arguments(boost::program_options::variables_map &vm);

  virtual void setup();

  virtual bool run();
};

#endif /* end of include guard: LENNARD_JONES_SYSTEM_H_ZDNDUQM8 */

/* vim: set ts=2 sw=2 sts=2 tw=0 expandtab :*/
