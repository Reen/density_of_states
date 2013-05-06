#include "ising_system.h"

// C++ Standard Library
#include <iomanip>

// Boost Bind
#include <boost/bind.hpp>

#include "mc_sampler.h"

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

/**
 * Public Methods
 */

/**
 * Constructor
 */
IsingSystem::IsingSystem() {}

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

void IsingSystem::setup(settings_t s) {
	SimulationSystem::setup(s);
	
	try {
		sampler       = boost::any_cast<size_t>(settings["sampler"]);
	} catch(const boost::bad_any_cast &e) {
		std::cerr << "Error while extracting settings in "
			<< __FILE__ << ":" << __LINE__ << "\nExecption: "
			<< e.what() << std::endl;
		throw e;
	}

	size_t error_acc_size = 9 * (size_t)log10(steps / error_check_f) + 1;
	error_acc.resize(error_acc_size);

	// set_size should have been called at this point
	// and thus n_bins should have a value
	settings["macrostates"] = n_bins;

	read_exact_dos();
	setup_output();
}

bool IsingSystem::run() {
	switch (sampler) {
	case 0:
		mc_loop<BoltzmannSampler<BoltzmannTableFunctor> >();
		break;
	case 1:
		mc_loop<WangLandauSampler>();
		break;
	/*case 2:
		mc_loop<QualityMeasureBSampler>();
		break;
	case 3:
		mc_loop<QualityMeasureASampler>();
		break;*/
	default:
		std::cerr << "Error: unknown sampler" << std::endl;
		return false;
	}

	out << "#\n#          time     mean_error      var_error   mean_error_s    var_error_s       qm_error\n";
	for (size_t i = 0; i < error_acc.size(); i++) {
		size_t time  = error_acc[i].step;
		double mean  = boost::accumulators::mean(error_acc[i].err1);
		double var   = boost::accumulators::variance(error_acc[i].err1);
		double mean2 = boost::accumulators::mean(error_acc[i].err2);
		double var2  = boost::accumulators::variance(error_acc[i].err2);
		//double qm_err   = error_acc[i].get<5>() / runs;
		out << std::setw(15) << std::right
			<< time
			<< std::setw(15) << std::right
			<< mean
			<< std::setw(15) << std::right
			<< var
			<< std::setw(15) << std::right
			<< mean2
			<< std::setw(15) << std::right
			<< var2
			//<< std::setw(15) << std::right
			//<< qm_err
			<< "\n";

	}
	return true;
}
