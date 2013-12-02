#include "simulation_rng.h"

#include <stdio.h>

size_t SimulationRNG::initialize() {
	if (!seed_set) {
		FILE* devran = fopen("/dev/urandom", "rb");
		size_t read = fread(&seed, sizeof(size_t), 1, devran);
		if (read != 1) {
			throw std::runtime_error("Could not read from /dev/urandom");
		}
		fclose(devran);
	}
	rng.seed(seed);
	seed_set = true;
	return seed;
}
