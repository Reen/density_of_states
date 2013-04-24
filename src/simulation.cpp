#include "simulation.h"

// Boost Program Options
#include <boost/program_options.hpp>

namespace po = boost::program_options;


void Simulation::parse_arguments(int argc, char *argv[]) {
	// parse argc / argv
	po::options_description all("Allowed options");
	po::options_description desc("General options");
	desc.add_options()
		("help", "produce help message")
		("steps,S",       po::value<size_t>()->default_value(1000000),  "Number of steps per simulation")
		("runs,R",        po::value<size_t>()->default_value(1000),      "Number of simulations")
		("temperature,T", po::value<double>()->default_value(2.0),          "Temperature")
		("tag",           po::value<std::string>(),                       "Additional tag to append to files")
		("seed",          po::value<size_t>(),                                "Random seed")
		("system",        po::value<std::string>(&system)->default_value("Ising"), "System to run")
		("flatness,f",    po::value<double>()->default_value(0.99),  "Flatness parameter of the Wang-Landau algorithm")
		;
	all.add(desc);
	for (simulation_systems_t::iterator iter = simulation_systems.begin();
					iter != simulation_systems.end();
					++iter) {
		all.add(iter->second->get_program_options());
	}
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, all), vm);
	po::notify(vm);

	if (vm.count("help")) {
		std::cout << "Ising Simulator "
			<< VERSION_MAJOR << "." << VERSION_MINOR
			<< "-" << g_GIT_SHA1 << std::endl;
		std::cout << all << std::endl;
		std::exit(EXIT_FAILURE);
	}

	if (vm.count("system")) {
		if (simulation_systems.count(vm["system"].as<std::string>()) != 1) {
			std::cerr << "Unknown System! Valid values are:\n";
			for (simulation_systems_t::iterator iter = simulation_systems.begin();
					iter != simulation_systems.end();
					++iter) {
				std::cerr << iter->first << "\n";
			}
			std::cerr << std::endl;
			std::exit(EXIT_FAILURE);
		}
	}

	settings["temperature"] = vm["temperature"].as<double>();
	settings["steps"]       = vm["steps"].as<size_t>();
	settings["runs"]        = vm["runs"].as<size_t>();
	settings["kB"]          = kB;
	settings["flatness"]    = vm["flatness"].as<double>();

	// safety check
	if (boost::any_cast<size_t>(settings["steps"]) % error_check_f != 0) {
		std::cerr << "Error: check-frequency and number of steps don't match" << std::endl;
		std::exit(EXIT_FAILURE);
	}

}

void Simulation::setup_variables() {

}

bool Simulation::run() {
	return simulation_systems[system]->run();
}

int Simulation::exec(int argc, char *argv[]) {
	try {
		parse_arguments(argc, argv);
		setup_variables();
	} catch(std::exception& e) {
		std::cerr << "error: " << e.what() << std::endl;
		return EXIT_FAILURE;
	} catch(...) {
		std::cerr << "Exception of unknown type!" << std::endl;
		return EXIT_FAILURE;
	}

	
	return (run() ? EXIT_SUCCESS : EXIT_FAILURE);
}

const double Simulation::kB = 1.0;

