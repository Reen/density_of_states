#ifndef SIMULATION_SYSTEM_H

#define SIMULATION_SYSTEM_H


// Boost Program Options
#include <boost/program_options.hpp>

#include "simulation_rng.h"
#include "typedefs.h"

class SimulationSystem {
protected:
	settings_t settings;
	size_t steps;
	size_t runs;
	size_t error_check_f;

	// transition counting matrix
	matrix_int_t Q;

public:
	virtual boost::program_options::options_description get_program_options() = 0;
	virtual bool run() = 0;
	virtual void setup(settings_t s) = 0;

	SimulationRNG rng;

};

#endif /* end of include guard: SIMULATION_SYSTEM_H */
