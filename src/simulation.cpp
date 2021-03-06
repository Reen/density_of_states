#include "simulation.h"

// Boost Program Options
#include <boost/program_options.hpp>

// for get_executable_path
#include<libgen.h>
#ifdef HAVE_DARWIN
#include <mach-o/dyld.h>
#endif

namespace po = boost::program_options;


void Simulation::parse_arguments(int argc, char *argv[]) {
	// rebuild the command line for documentation
	std::string cmdline = "";
	for (int i = 0; i < argc; ++i)
	{
		cmdline += " ";
		cmdline += argv[i];
	}

	// parse argc / argv
	po::options_description all("Allowed options");
	po::options_description desc("General options");
	desc.add_options()
		(
			"help",
			"produce help message"
		) (
			"steps,S",
			po::value<size_t>()->default_value(1000000),
			"Number of steps per simulation"
		) (
			"runs,R",
			po::value<size_t>()->default_value(1000),
			"Number of simulations"
		) (
			"temperature,T",
			po::value<double>()->default_value(2.0),
			"Temperature"
		) (
			"tag",
			po::value<std::string>()->default_value(""),
			"Additional tag to append to files"
		) (
			"seed",
			po::value<size_t>(),
			"Random seed"
		) (
			"system",
			po::value<std::string>(&system)->default_value("Ising"),
			"System to run, choose:\nIsing\nToyDos\nLJones"
		) (
			"sampler",
			po::value<size_t>()->default_value(1),
			"Sampler to use:\n0 – Boltzmann, 1 – Wang-Landau, 2 – Q-Matrix B, 3 – Q-Matrix A, 4 – TransitionMatrix, 5 – 1/t Algorithm"
		) (
			"flatness,f",
			po::value<double>()->default_value(0.99),
			"Flatness parameter of the Wang-Landau algorithm"
		) (
			"error-check-freq",
			po::value<size_t>()->default_value(100),
			"Frequency of error checking"
		) (
			"schedule",
			po::value<size_t>()->default_value(0),
			"Temperature Schedule to use for Sampler 0:\n0 – none, 1 – Linear, 2 – Sine,\n3 – LogSine, 4 – Exp, 5 – Pow,\n6 – Log"
		) (
			"schedule-Tmin",
			po::value<double>()->default_value(0.1),
			"Minimum temperature for schedules."
		) (
			"schedule-A",
			po::value<double>()->default_value(5.0),
			"Number of oscillations in Sine and LogSine schedule."
		) (
			"schedule-alpha",
			po::value<double>()->default_value(1.0-1e-8),
			"Exponential schedule: term alpha in alpha^t"
		) (
			"schedule-B",
			po::value<double>()->default_value(10),
			"Power law schedule: T0 * B * t**gamma"
		) (
			"schedule-gamma",
			po::value<double>()->default_value(-0.39),
			"Power law schedule: T0 * B * t**gamma"
		) (
			"one-over-t-c",
			po::value<double>()->default_value(1.0),
			"Parameter c of 1/t alg., i.e. c/t"
		) (
			"one-over-t-s",
			po::value<size_t>()->default_value(1500),
			"Parameter s of 1/t alg., i.e. s * macrostates"
		);
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
		} else {
			simulation_systems[vm["system"].as<std::string>()]->parse_arguments(vm);
		}
	}

	settings["temperature"]   = vm["temperature"].as<double>();
	settings["steps"]         = vm["steps"].as<size_t>();
	settings["runs"]          = vm["runs"].as<size_t>();
	settings["kB"]            = kB;
	settings["flatness"]      = vm["flatness"].as<double>();
	settings["error_check_f"] = vm["error-check-freq"].as<size_t>();
	settings["tag"]           = vm["tag"].as<std::string>();
	settings["sampler"]       = vm["sampler"].as<size_t>();
	settings["cmdline"]       = cmdline;
	settings["schedule"]      = vm["schedule"].as<size_t>();
	settings["schedule-Tmin"] = vm["schedule-Tmin"].as<double>();
	settings["schedule-A"]    = vm["schedule-A"].as<double>();
	settings["schedule-alpha"]= vm["schedule-alpha"].as<double>();
	settings["schedule-B"]    = vm["schedule-B"].as<double>();
	settings["schedule-gamma"]= vm["schedule-gamma"].as<double>();
	settings["one-over-t-c"]  = vm["one-over-t-c"].as<double>();
	settings["one-over-t-s"]  = vm["one-over-t-s"].as<size_t>();

}

std::string get_executable_path() {
	// determin path this executable lives in
	// to add more OS see:
	// http://stackoverflow.com/questions/1023306/finding-current-executables-path-without-proc-self-exe
	// http://stackoverflow.com/questions/799679/programatically-retrieving-the-absolute-path-of-an-os-x-command-line-app
	char dir[1024] = ".";
#ifdef HAVE_LINUX
	ssize_t len = readlink("/proc/self/exe", dir, sizeof(dir)-1);
	if (len != -1) {
		dir[len] = '\0';
		dirname(dir); // dirname modifies argument
	} else {
		throw std::runtime_error("An error occured while determining the executable path");
	}
#endif
#ifdef HAVE_DARWIN
	char buf[1024];
	uint32_t size = sizeof(buf);
	if (_NSGetExecutablePath(buf, &size) == 0) {
		realpath(buf, dir);
		// dirname does not modify argument on OSX but instead returns a pointer
		// to internal memory. See man 3 dirname
		strcpy(dir, dirname(dir));
	} else {
		throw std::runtime_error("An error occured while determining the executable path");
	}
#endif
#if VERBOSE==1
	std::cout << "Path to executable: " << dir << std::endl;
#endif
	return std::string(dir);
}

void Simulation::setup_variables() {
	settings["executable_path"] = get_executable_path();

	simulation_systems[system]->setup();
}

bool Simulation::run() {
	return simulation_systems[system]->run();
}

int Simulation::exec(int argc, char *argv[]) {

	parse_arguments(argc, argv);
	setup_variables();

	return (run() ? EXIT_SUCCESS : EXIT_FAILURE);
}

const double Simulation::kB = 1.0;

