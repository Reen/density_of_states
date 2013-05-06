#ifndef ISING_SYSTEM_H

#define ISING_SYSTEM_H

// C++ Standard Library
#include <vector>

// Boost Accumulator
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/variance.hpp>

// Boost Assert
#include <boost/assert.hpp>

// Boost Format
#include <boost/format.hpp>

// Boost Program Options
#include <boost/program_options.hpp>

// Boost Random
#include <boost/random/uniform_01.hpp>


#include "simulation_system.h"
#include "q_matrix_tools.h"


using namespace rhab;

struct StepStatistics {
	size_t step;
	boost::accumulators::accumulator_set<
		double,
		boost::accumulators::stats<
			boost::accumulators::tag::min,
			boost::accumulators::tag::max,
			boost::accumulators::tag::mean,
			boost::accumulators::tag::variance,
			boost::accumulators::tag::median
				> > err1;
	boost::accumulators::accumulator_set<
		double,
		boost::accumulators::stats<
			boost::accumulators::tag::min,
			boost::accumulators::tag::max,
			boost::accumulators::tag::mean,
			boost::accumulators::tag::variance,
			boost::accumulators::tag::median
				> > err2;
};

class IsingSystem : public SimulationSystem {
private:
	typedef boost::numeric::ublas::matrix<signed char> storage_t;
	typedef std::vector< StepStatistics > error_acc_t;

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

	// error accumulator array
	error_acc_t error_acc;

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
			std::cout << run << std::endl;
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
					// calculate_statistics()
					error_acc[index].step = step;
					Qd = normalize_q(Q);

					double err = calculate_error_q(dos_exact_norm, Qd);
					error_acc[index].err1(err);
					if (sampler.has_own_statistics()) {
						err = sampler.calculate_error(dos_exact_norm);
						error_acc[index].err2(err);
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


	void setup_output() {
		std::vector<std::string> sampler_string;
		// initialize list of available samplers
		sampler_string.push_back("BM");
		sampler_string.push_back("WL");
		sampler_string.push_back("QB");
		sampler_string.push_back("QA");

		std::string format;
		switch (sampler) {
			case 1:
				format = "dos_%1%_%2%S_%3%R_%4%M_%6$0.2ff%7%%8%.out";
				break;
			case 2:
			case 3:
				format = "dos_%1%_%2%S_%3%R_%4%M%7%%8%.out";
				break;
			case 0:
			default:
				format = "dos_%1%_%2%S_%3%R_%4%M_%5$0.2fT%7%%8%.out";
		}
		std::string buf = str( boost::format(format)
				% sampler_string[sampler]
				% steps
				% runs
				% n_bins
				% boost::any_cast<double>(settings["temperature"])
				% boost::any_cast<double>(settings["flatness"])
				% (tag.size() > 0 ? "_" : "")
				% tag
			);
		out.open( buf.c_str() );
		if (!out.good()) {
			throw std::runtime_error("Error: could not open output file");
		}
	}

public:
	IsingSystem();

	virtual boost::program_options::options_description get_program_options();

	virtual void setup(settings_t s);

	virtual bool run();
};

#endif /* end of include guard: ISING_SYSTEM_H */
