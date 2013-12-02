#ifndef SIMULATION_RNG_H

#define SIMULATION_RNG_H
//
// Boost Random
#include <boost/random/mersenne_twister.hpp>

struct SimulationRNG {
	// random number generator
	boost::mt19937 rng;
	size_t seed;
	bool seed_set;

	SimulationRNG()
		: seed_set(false) {}
	size_t initialize();
};

#endif /* end of include guard: SIMULATION_RNG_H */
