#ifndef SIMULATION_H

#define SIMULATION_H

// C++ Standard Library
#include <map>
#include <string>
#include <vector>

// Boost Assign
#include <boost/assign/list_inserter.hpp>
#include <boost/assign/std/vector.hpp>

// Boost Smart Pointer
#include <boost/shared_ptr.hpp>

// Project Includes
#include "config.h"
#include "GitSHA1.h"
#include "simulation_system.h"
#include "ising_system.h"
#include "toydos_system.h"



class Simulation {
private:
	typedef std::map< std::string, boost::shared_ptr<SimulationSystem> > simulation_systems_t;
	// variables set in constructor
	std::vector< std::string > sampler;
	simulation_systems_t simulation_systems;

	// variables set by parse_arguments
	settings_t settings;
	std::string tag;
	std::string system;

	/**
	 * Methods
	 */
	void parse_arguments(int argc, char *argv[]);

	void setup_variables();

	bool run();

public:
	static const double kB;

	Simulation ()
	{
		using namespace boost::assign;
		insert(simulation_systems)
			("Ising",  boost::shared_ptr<SimulationSystem>(new IsingSystem(settings)))
			("ToyDos", boost::shared_ptr<SimulationSystem>(new ToyDosSystem(settings)))
			;
	}

	int exec(int argc, char *argv[]);
};



#endif /* end of include guard: SIMULATION_H */
