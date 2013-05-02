#ifndef ISING_SYSTEM_H

#define ISING_SYSTEM_H

// C++ Standard Library
#include <vector>

// Boost Bind
#include <boost/bind.hpp>

// Boost Program Options
#include <boost/program_options.hpp>

// Boost Random
#include <boost/random/uniform_01.hpp>


#include "simulation_system.h"
#include "mc_sampler.h"

namespace po = boost::program_options;


class IsingSystem : public SimulationSystem {
private:
	typedef boost::numeric::ublas::matrix<signed char> storage_t;

	size_t size;
	storage_t lattice;
	size_t sampler;
	double J;
	boost::uniform_01<> dist01;

	// size dependent constants:
	int e_min;
	int e_max;
	int n_bins;

	void set_size(size_t L) {
		size = L;
		lattice.resize(L, L);
		e_min = -2 * L * L;
		e_max =  2 * L * L;
		n_bins = e_max - e_min;
		Q.resize(n_bins, n_bins);
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

				// update Q matrix
				Q(energy-e_min, energy-e_min+dE)++;

				if (sampler(energy, energy+dE, energy, energy+dE)) {
					lattice(i, j) *= -1;
					energy += dE;
				}
				if (step % error_check_f == 0) {
					// calculate_statistics()
				}
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
			("ising-sampler",     po::value<size_t>(&sampler)->default_value(1),      "Sampler to use:\n0 – Boltzmann, 1 – Wang-Landau, 2 – Q-Matrix")
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
			// safety check
			if (steps % error_check_f != 0) {
				std::cerr << "Error: check-frequency and number of steps don't match" << std::endl;
				std::exit(EXIT_FAILURE);
			}
		} catch(const boost::bad_any_cast &e) {
			std::cerr << "bla: " << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}


	}

	virtual bool run() {
		switch (sampler) {
		case 0:
			mc_loop<BoltzmannSampler>();
			break;
		/*case 1:
			mc_loop<WangLandauSampler>();
			break;
		case 2:
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
