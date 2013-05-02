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

// Boost Bind
#include <boost/bind.hpp>

// Boost Program Options
#include <boost/program_options.hpp>

// Boost Random
#include <boost/random/uniform_01.hpp>


#include "simulation_system.h"
#include "mc_sampler.h"
#include "q_matrix_tools.h"

namespace po = boost::program_options;

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
	int n_bins;

	void set_size(size_t L) {
		size = L;
		lattice.resize(L, L);
		e_min = -2 * L * L;
		e_max =  2 * L * L;
		n_bins = e_max - e_min; // evtl durch 4 teilen
		Q.resize(n_bins, n_bins);
		Qd.resize(n_bins, n_bins);
	}

	int reset() {
		int magnetization = 0;
		// set lattice to random state
		for (size_t i = 0; i < lattice.data().size(); ++i) {
			lattice.data()[i] = (dist01(rng.rng) > 0.5 ? 1 : -1);
			magnetization += lattice.data()[i];
		}
		return magnetization;
	}

	int calculate_energy() {
		int energy = 0;
		for (size_t i = 0; i < size; ++i) {
			for (size_t j = 0; j < size; ++j) {
				energy += -J * (lattice(i, j) * (
							lattice((i-1+size)%size,  j              ) +
							lattice( i             , (j-1+size)%size )
							));
			}
		}
		return energy;
	}

	template<class Sampler>
	void mc_loop() {
		
		for (size_t run = 0; run < runs; run++) {
			// variables for error / statistics calculation
			size_t error_check_freq = error_check_f;
			size_t index = 0;
			// reset Q matrix
			Q *= 0;
			int magnetization = reset();
			int energy        = calculate_energy();
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

				// update Q matrix
				Q(Ei, Ej)++;

				if (sampler(energy, energy+dE, Ei, Ej)) {
					lattice(i, j) *= -1;
					energy += dE;
				}
				if (step % error_check_f == 0) {
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
					if (step % (10*error_check_f) == 0) {
						error_check_f *= 10;
					}
				}
				sampler.check(step, run);
			}
		}
	}

public:
	IsingSystem()
	{}

	virtual boost::program_options::options_description get_program_options() {
		boost::program_options::options_description desc("Ising System Options");
		desc.add_options()
			("ising-size",        po::value<size_t>(&size)->default_value(50)->notifier(boost::bind(&IsingSystem::set_size, this, _1)), "Number of Spins in one dimension.\nGrid will be size*size")
			("ising-sampler",     po::value<size_t>(&sampler)->default_value(1),      "Sampler to use:\n0 – Boltzmann, 1 – Wang-Landau, 2 – Q-Matrix B, 3 – Q-Matrix A")
			("ising-interaction", po::value<double>(&J)->default_value(1),      "Interaction J\nJ>0 - ferromagnetic\nJ<0 - antiferromagnetic\nJ=0 - noninteracting")
			;
		return desc;
	}

	virtual void setup(settings_t s) {
		settings = s;

		set_size(size);

		if (settings.count("seed")) {
			rng.seed = boost::any_cast<size_t>(settings["seed"]);
			rng.seed_set = true;
		} else {
			rng.initialize();
		}

		try {
			runs          = boost::any_cast<size_t>(settings["runs"]);
			steps         = boost::any_cast<size_t>(settings["steps"]);
			error_check_f = boost::any_cast<size_t>(settings["error_check_f"]);
		} catch(const boost::bad_any_cast &e) {
			std::cerr << "bla: " << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}

		// safety check
		if (steps % error_check_f != 0) {
			std::cerr << "Error: check-frequency and number of steps don't match" << std::endl;
			std::exit(EXIT_FAILURE);
		}
		size_t error_acc_size = 9 * (size_t)log10(steps / error_check_f) + 1;
		error_acc.resize(error_acc_size);


	}

	virtual bool run() {
		switch (sampler) {
		case 0:
			mc_loop<BoltzmannSampler>();
			break;
		case 1:
			mc_loop<WangLandauSampler>();
			break;
		/*case 2:
			mc_loop<QualityMeasureBSampler>();
			break;
		case 3:
			mc_loop<QualityMeasureASampler>();
			break;*/
		default:
			std::cerr << "Error: unknown sampler" << std::endl;
			return false;
		}
		return true;
	}
};

#endif /* end of include guard: ISING_SYSTEM_H */
