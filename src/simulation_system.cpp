#include "simulation_system.h"

void SimulationSystem::setup(settings_t s) {
	settings = s;

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
		tag           = boost::any_cast<std::string>(settings["tag"]);
	} catch(const boost::bad_any_cast &e) {
		std::cerr << "Error while extracting settings in "
			<< __FILE__ << ":" << __LINE__ << "\nExecption: "
			<< e.what() << std::endl;
		throw e;
	}

	// safety check
	if (steps % error_check_f != 0) {
		std::cerr << "Error: check-frequency and number of steps don't match" << std::endl;
		std::exit(EXIT_FAILURE);
	}
}
