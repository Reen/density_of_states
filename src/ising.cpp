// C++ Standard Library
#include <map>
#include <string>
#include <vector>

// Boost Assign
#include <boost/assign/list_inserter.hpp>
#include <boost/assign/std/vector.hpp>

// Boost Bind
#include <boost/bind.hpp>

// Boost Program Options
#include <boost/program_options.hpp>

// Boost Random
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_01.hpp>

// Boost Smart Pointer
#include <boost/shared_ptr.hpp>


// Project Includes
#include "config.h"
#include "GitSHA1.h"

namespace po = boost::program_options;

class SimulationSystem {
public:
	virtual boost::program_options::options_description get_program_options() = 0;
};

class IsingSystem : public SimulationSystem {
	typedef std::vector<bool> storage_t;
private:
	size_t size;
	storage_t spins;

	void set_size(size_t L) {
		size = L;
		spins.resize(L*L);
	}

public:
	IsingSystem()
	{}

	virtual boost::program_options::options_description get_program_options() {
		boost::program_options::options_description desc("Ising System Options");
		desc.add_options()
			("ising-size", po::value<size_t>(&size)->default_value(50)->notifier(boost::bind(&IsingSystem::set_size, this, _1)), "Number of Spins in one dimension.\nGrid will be size*size")
			;
		return desc;
	}
};

class Simulation {
private:
	typedef std::map< std::string, boost::shared_ptr<SimulationSystem> > simulation_systems_t;
	// variables set in constructor
	std::vector< std::string > sampler;
	simulation_systems_t simulation_systems;
	size_t error_check_f;

	// variables set by parse_arguments
	size_t steps;
	size_t runs;
	std::string tag;
	std::string system;

	// random number generator
	boost::mt19937 rng;
	size_t seed;
	bool seed_set;

	void parse_arguments(int argc, char *argv[]) {
		// parse argc / argv
		po::options_description all("Allowed options");
		po::options_description desc("General options");
		desc.add_options()
			("help", "produce help message")
			("steps,S",       po::value<size_t>(&steps)->default_value(1000000),  "Number of steps per simulation")
			("runs,R",        po::value<size_t>(&runs)->default_value(1000),      "Number of simulations")
			//("temperature,T", po::value<double>(&T)->default_value(2.0),          "Temperature")
			("tag",           po::value<std::string>(&tag),                       "Additional tag to append to files")
			("seed",          po::value<size_t>(),                                "Random seed")
			("system",        po::value<std::string>(&system)->default_value("Ising"), "System to run")
			//("sampler",       po::value<size_t>(&sampler)->default_value(1),      "Sampler to use: 0 – Boltzmann, 1 – Wang-Landau, 2 – Q-Matrix")
			//("flatness,f",    po::value<double>(&flatness)->default_value(0.99),  "Flatness parameter of the Wang-Landau algorithm")
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

		if (vm.count("seed")) {
			seed = vm["seed"].as<size_t>();
			seed_set = true;
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
	}

	void initialize_rng()
	{
		if (!seed_set) {
			FILE* devran = fopen("/dev/urandom", "rb");
			size_t read = fread(&seed, sizeof(size_t), 1, devran);
			if (read != 1) {
				throw std::runtime_error("Could not read from /dev/urandom");
			}
			fclose(devran);
		}
		rng.seed(seed);
	}

	void safety_check()
	{
	}

	void setup_variables()
	{
	}

	bool run() {

		return true;
	}

public:
	static const double kB;

	Simulation ()
		: seed_set(false)
	{
		using namespace boost::assign;
		insert(simulation_systems)
			("Ising", boost::shared_ptr<SimulationSystem>(new IsingSystem()))
			;
	}

	int exec(int argc, char *argv[]) {
		try {
			parse_arguments(argc, argv);
			initialize_rng();
			setup_variables();
			safety_check();
		} catch(std::exception& e) {
			std::cerr << "error: " << e.what() << std::endl;
			return EXIT_FAILURE;
		} catch(...) {
			std::cerr << "Exception of unknown type!" << std::endl;
			return EXIT_FAILURE;
		}

		// safety check
		if (steps % error_check_f != 0) {
			std::cerr << "Error: check-frequency and number of steps don't match" << std::endl;
			return EXIT_FAILURE;
		}
		
		return (run() ? EXIT_SUCCESS : EXIT_FAILURE);
	}
};

const double Simulation::kB = 1.0;

int main(int argc, char *argv[])
{
	Simulation s;
	return s.exec(argc, argv);
}
