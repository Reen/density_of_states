#ifndef ISING_SYSTEM_H

#define ISING_SYSTEM_H

// C++ Standard Library
#include <vector>

// Boost Program Options
#include <boost/program_options.hpp>

// Boost Random
#include <boost/random/uniform_01.hpp>


#include "simulation_system.h"
#include "q_matrix_tools.h"


//using namespace rhab;


class IsingSystem : public SimulationSystem {
private:
	typedef boost::numeric::ublas::matrix<signed char> storage_t;

	// size L of one dimension of the lattice
	size_t size;

	// the lattice itself
	storage_t lattice;

	/**
	 * Settings
	 */
	// the sampler chosen
	size_t sampler;

	// the interaction constant J
	double J;


	/**
	 * Misc
	 */
	boost::uniform_01<> dist01;

	// exact dos
	vector_double_t dos_exact_norm;

	// size dependent constants:
	int e_min;
	int e_max;
	size_t n_bins;

	void set_size(size_t L);

	int reset();

	int calculate_energy();

	template<class Sampler>
	void mc_loop() {
		for (size_t run = 0; run < runs; run++) {
			std::cout << "run: " << run << std::endl;
			// variables for error / statistics calculation
			size_t error_check_freq = error_check_f;
			size_t index = 0;
			// reset Q matrix
			Q *= 0;
			int magnetization(0), energy(1);
			while (energy > 0) {
				magnetization = reset();
				energy        = calculate_energy();
				(void)magnetization;
			}
			int i, j;

			Sampler sampler(rng.rng, Q, settings);
			for (size_t step = 1; step <= steps; step++) {
				// single spin flip:
				i = (int)(dist01(rng.rng) * size);
				j = (int)(dist01(rng.rng) * size);

				// calculate change in energy
				int dE = 2 * J * lattice(i,j) * (
						lattice((i+1)%size     ,  j             ) +
						lattice((i-1+size)%size,  j             ) +
						lattice( i             , (j+1)%size     ) +
						lattice( i             , (j-1+size)%size)
						);

				BOOST_ASSERT((energy-e_min)%4 == 0);
				BOOST_ASSERT((energy-e_min+dE)%4 == 0);
				int Ei = (energy-e_min) / 4;
				int Ej = (energy-e_min+dE) / 4;
				BOOST_ASSERT(Ei >= 0);
				BOOST_ASSERT(Ej >= 0);

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
					lattice(i, j) *= -1;
					energy += dE;
				}
				if (step % error_check_freq == 0) {
					error_acc[index].step = step;

					double err_lq, err_gth, err_power;
					bool   suc_lq, suc_gth, suc_power;
					boost::tie(err_lq, err_gth, err_power, suc_lq, suc_gth, suc_power) = rhab::calculate_error_q(dos_exact_norm, Q, Qd);
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
						double err = sampler.calculate_error(dos_exact_norm);
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

	void read_exact_dos();

	void setup_output();

public:
	IsingSystem(settings_t & s);

	virtual boost::program_options::options_description get_program_options();
	virtual void parse_arguments(boost::program_options::variables_map &vm);

	virtual void setup();

	virtual bool run();
};

#endif /* end of include guard: ISING_SYSTEM_H */
