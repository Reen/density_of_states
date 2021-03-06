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

	void set_size(size_t L);

	int reset();

	int calculate_energy();

	template<class Sampler>
	void mc_loop() {
		size_t index2 = 1;

		vector_double_t dos_lsq(Q.size1());
		vector_double_t dos_gth(Q.size1());
		vector_double_t dos_pow(Q.size1());

		size_t run_from;
		size_t run_to;
		assert(runs % world_size == 0);
		run_from = (runs/world_size) *  world_rank;
		run_to   = (runs/world_size) * (world_rank + 1);

#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif
		for (size_t run = run_from; run < run_to; run++) {
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
					error_acc.set_step(index, step);

					double err_lsq, err_gth, err_pow;
					bool   suc_lsq, suc_gth, suc_pow;
					double err_qm;
					boost::tie(
							err_lsq, err_gth, err_pow,
							suc_lsq, suc_gth, suc_pow, err_qm) =
						rhab::calculate_error_q(
								dos_exact_norm, matrix_double_t(), Q, Qd, error_matrices,
								dos_lsq, dos_gth, dos_pow, index);
					if (suc_lsq) {
						error_acc.push(0, index, err_lsq);
					}
					if (suc_gth) {
						error_acc.push(1, index, err_gth);
					}
					if (suc_pow) {
						error_acc.push(2, index, err_pow);
					}

					if (sampler.has_own_statistics()) {
						double err = sampler.calculate_error(dos_exact_norm, error_matrices, index);
						error_acc.push(3, index, err);
					}

					// we only capture the parameter for the first run
					if (sampler.has_parameter() && run == 0) {
						double param(0.0);
						sampler.get_parameter(param);
						error_acc.set_parameter(index, param);
					}

					index++;
					if (step % (10*error_check_freq) == 0) {
						error_check_freq *= 10;
					}
				}
				sampler.check(step, run);
			}

			if (world_size == 1 && run+1 == index2) {
				std::ostringstream add;
				add << "# last Q/Qd matrix:\n# " << Q << "\n# " << Qd << std::endl;
				write_output(run+1, add.str());
				index2 *= 10;
			}
		}
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		error_acc.pull();
		combine_err_matrix(error_matrices);
		if (world_rank == 0) {
			std::ostringstream add;
			add << "# last Q/Qd matrix:\n# " << Q << "\n# " << Qd << std::endl;
			write_output(runs, add.str());
		}
#endif
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

