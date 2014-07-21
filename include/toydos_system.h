#ifndef TOYDOS_SYSTEM_H_APVVKHXZ
#define TOYDOS_SYSTEM_H_APVVKHXZ

// Boost Random
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int.hpp>

#include "simulation_system.h"
#include "q_matrix_tools.h"

class ToyDosSystem : public SimulationSystem {
private:
	typedef boost::numeric::ublas::matrix<signed char> storage_t;

	/**
	 * Settings
	 */
	// the sampler chosen
	size_t sampler;
	size_t macrostates;
	size_t connections;

	/**
	 * Toy System Properties
	 */
	int gengraph_seed;
	bool gengraph_seed_set;
	void set_gengraph_seed(int);
	vector_int_t get_exact_dos(size_t macro_states);
	vector_int_t get_energy_map(const vector_int_t &dos_exact, const size_t &num_config);
	bool has_regular_graph(size_t n, size_t r);


	// exact density of states
	vector_int_t dos_exact;

	// number of configurations
	size_t num_config;

	// exact density of states normalized
	vector_double_t dos_exact_norm;

	// map configuration to energy
	vector_int_t config_to_energy;

	// connectivity of states
	matrix_int_t mt;

	// microstate transition matrix
	matrix_double_t microstate_tm;

	// path to graph utility
	std::string graph_bin;

	// lumping matrix
	matrix_int_t lumper;

	// exact Q matrix
	matrix_double_t exact_q;

	void setup_variables();
	void safety_check();

	/**
	 * Misc
	 */
	boost::uniform_01<> dist01;

	// size dependent constants:
	int e_min;
	int e_max;

	void set_size(size_t L);

	template<class Sampler>
	void mc_loop() {
		boost::uniform_int<> select_pos(0,connections-1);
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
			// start at random position
			size_t state = select_pos(rng.rng);
			// reset Q matrix
			Q *= 0;
			Sampler sampler(rng.rng, Q, settings);
			for (size_t step = 1; step <= steps; step++) {
				size_t new_state = mt(state, select_pos(rng.rng));
				int Eold = config_to_energy[state];
				int Enew = config_to_energy[new_state];
				Q(Eold, Enew)++;
				assert(new_state != state);
				if (sampler(Enew-Eold, Eold, Enew)) {
					state = new_state;
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
								dos_exact_norm, exact_q, Q, Qd, error_matrices,
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

					error_acc.push(4, index, err_qm);

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

	void setup_output();

public:
	ToyDosSystem(settings_t &s);

	virtual boost::program_options::options_description get_program_options();

	virtual void parse_arguments(boost::program_options::variables_map &vm);

	virtual void setup();

	virtual bool run();

};
#endif /* end of include guard: TOYDOS_SYSTEM_H_APVVKHXZ */

