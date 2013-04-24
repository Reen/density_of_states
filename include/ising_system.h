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
	typedef std::vector<bool> storage_t;
private:
	size_t size;
	storage_t spins;
	size_t sampler;
	double J;
	boost::uniform_01<> dist01;

	void set_size(size_t L) {
		size = L;
		spins.resize(L*L);
	}

	void reset() {
		// set spins to random state
		for (size_t i = 0; i < spins.size(); ++i) {
			spins[i] = (dist01(rng.rng) > 0.5);
		}
	}

	void energy() {
		
	}

	template<class Sampler>
	void mc_loop() {
		
		for (size_t run = 0; run < runs; run++) {
			reset();

			Sampler sampler(rng.rng, Q, settings);
			for (size_t step = 1; step <= steps; step++) {
				
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
			runs    = boost::any_cast<size_t>(settings["runs"]);
			steps   = boost::any_cast<size_t>(settings["steps"]);
		} catch(const boost::bad_any_cast &e) {
			std::cerr << "bla: " << e.what() << std::endl;
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
