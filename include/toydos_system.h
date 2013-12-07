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
		for (size_t run = 0; run < runs; run++) {
			error_check_f = 100;
			// start at random position
			size_t state = select_pos(rng.rng);
			Q *= 0;
			//std::cout << "Q: " << Q << std::endl;
			size_t step = 1;
			size_t index = 0;
			Sampler sampler(rng.rng, Q, settings);
			for (; step <= steps; step++) {
				size_t new_state = mt(state, select_pos(rng.rng));
				int Eold = config_to_energy[state];
				int Enew = config_to_energy[new_state];
				Q(Eold, Enew)++;
				assert(new_state != state);
				if (sampler(Enew-Eold, Eold, Enew)) {
					state = new_state;
				}
				if (step > 0 && step % error_check_f == 0) {
					error_acc[index].step = step;

					double err_lq, err_gth, err_power;
					bool   suc_lq, suc_gth, suc_power;
					boost::tie(err_lq, err_gth, err_power, suc_lq, suc_gth, suc_power) = rhab::calculate_error_q(dos_exact_norm, Q, Qd, error_matrices, index);
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

					// @todo
					//error_acc[index].get<5>() += rhab::calculate_error_q_matrix(exact_q, Qd);

					index++;
					if (step % (10*error_check_f) == 0) {
						error_check_f *= 10;
					}
				}
				sampler.check(step, run);
			}
			//std::cout << Q << std::endl;
			//std::cout << Qd << std::endl;
			//vector_double_t dos = calculate_dos_gth(Qd);
			//final_errors.push_back(calculate_error(dos_exact_norm, dos));
			//std::cout << run << " " << dos << " " << calculate_error(dos_exact_norm, dos) << " " << sampler.calculate_error(dos_exact_norm) << std::endl;
		}
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

