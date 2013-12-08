#include "ising_system.h"
#include "mc_sampler.h"

// C++ Standard Library
#include <iomanip>

// Boost Bind
#include <boost/bind.hpp>

// Boost Format
#include <boost/format.hpp>

namespace po = boost::program_options;


/**
 * Private Methods:
 */

void IsingSystem::set_size(size_t L) {
	size = L;
	lattice.resize(L, L);
	e_min = -2 * L * L;
	//e_max =  2 * L * L;
	e_max =  0;
	n_bins = (e_max - e_min)/4 + 1;
	Q.resize(n_bins, n_bins);
	Qd.resize(n_bins, n_bins);

	std::fill(Q.data().begin(), Q.data().end(), 0);
	std::fill(Qd.data().begin(), Qd.data().end(), 0.0);
}

int IsingSystem::reset() {
	int magnetization = 0;
	// set lattice to random state
	for (size_t i = 0; i < lattice.data().size(); ++i) {
		lattice.data()[i] = (dist01(rng.rng) > 0.5 ? 1 : -1);
		magnetization += lattice.data()[i];
	}
	return magnetization;
}

int IsingSystem::calculate_energy() {
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

void IsingSystem::read_exact_dos() {
	std::string format = "%1%/data/ising/exact_%2%.dat";
	std::string fn = str( boost::format(format)
			% boost::any_cast<std::string>(settings["executable_path"])
			% size
			);
	std::ifstream in(fn.c_str());

	dos_exact_norm.resize(n_bins);
	dos_exact_norm *= 0;
	for (size_t i = 0; i < n_bins; i++) {
		int energy;
		double dos;
		in >> energy >> dos;
		dos_exact_norm[i] = dos;
	}
}

void IsingSystem::setup_output() {
	std::vector<std::string> sampler_string;
	// initialize list of available samplers
	sampler_string.push_back("BM");
	sampler_string.push_back("WL");
	sampler_string.push_back("QB");
	sampler_string.push_back("QA");
	sampler_string.push_back("TM");

	std::string format;
	switch (sampler) {
		case 1:
			format = "ising_%1%_%2%S_%3%R_%4%M_%6$0.2ff%7%%8%.out";
			break;
		case 2:
		case 3:
		case 4:
			format = "ising_%1%_%2%S_%3%R_%4%M%7%%8%.out";
			break;
		case 0:
			format = "ising_%1%_%2%S_%3%R_%4%M_%5$0.2fT%7%%8%.out";
			break;
		default:
			throw std::runtime_error("Unknown smapler");
	}
	std::string buf = str( boost::format(format)
			% sampler_string[sampler]
			% steps
			% runs
			% n_bins
			% boost::any_cast<double>(settings["temperature"])
			% boost::any_cast<double>(settings["flatness"])
			% (tag.size() > 0 ? "_" : "")
			% tag
		);

	open_output_files(buf);

	write_header();

	out << "#";
	out << "\n# sampler:     " << sampler << " - " << sampler_string[sampler];
	out << "\n# steps:       " << steps;
	out << "\n# runs:        " << runs;
	out << "\n# n_bins:      " << n_bins;
	out << "\n# temperature: " << boost::any_cast<double>(settings["temperature"]);
	out << "\n# flatness:    " << boost::any_cast<double>(settings["flatness"]);
	out << "\n#" << std::endl;
}

/**
 * Public Methods
 */

/**
 * Constructor
 */
IsingSystem::IsingSystem(settings_t &s) : SimulationSystem(s) {}

po::options_description IsingSystem::get_program_options() {
	po::options_description desc("Ising System Options");
	desc.add_options()
		(
			"ising-size",
			po::value<size_t>(&size)->default_value(14)->notifier(boost::bind(&IsingSystem::set_size, this, _1)),
			"Number of Spins in one dimension.\nGrid will be size*size"
		) (
			"ising-interaction",
			po::value<double>(&J)->default_value(1),
			"Interaction J\nJ>0 - ferromagnetic\nJ<0 - antiferromagnetic\nJ=0 - noninteracting"
		);
	return desc;
}

void IsingSystem::parse_arguments(boost::program_options::variables_map &) {}

void IsingSystem::setup() {
	SimulationSystem::setup();

	try {
		sampler       = boost::any_cast<size_t>(settings["sampler"]);
	} catch(const boost::bad_any_cast &e) {
		std::cerr << "Error while extracting settings in "
			<< __FILE__ << ":" << __LINE__ << "\nExecption: "
			<< e.what() << std::endl;
		throw e;
	}

	// set_size should have been called at this point
	// and thus n_bins should have a value
	settings["macrostates"] = n_bins;

	read_exact_dos();
	setup_output();
}

bool IsingSystem::run() {
	switch (sampler) {
	case 0:
		mc_loop<BoltzmannSampler>();
		break;
	case 1:
		mc_loop<WangLandauSampler>();
		break;
	case 2:
		mc_loop<QualityMeasureBSampler>();
		break;
	case 3:
		mc_loop<QualityMeasureASampler>();
		break;
	case 4:
		mc_loop<TransitionMatrixSampler>();
		break;
	default:
		std::cerr << "Error: unknown sampler" << std::endl;
		return false;
	}

	write_output();

	out << "# last Q/Qd matrix:\n# " << Q << "\n# " << Qd << std::endl;

	return true;
}
