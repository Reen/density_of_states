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
		default:
			format = "ising_%1%_%2%S_%3%R_%4%M_%5$0.2fT%7%%8%.out";
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
	out.open( buf.c_str() );
	if (!out.good()) {
		throw std::runtime_error("Error: could not open output file");
	}
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

void IsingSystem::parse_arguments(boost::program_options::variables_map &vm) {}

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
	case 2:
		mc_loop<QualityMeasureBSampler>();
		break;
	case 3:
		mc_loop<QualityMeasureASampler>();
		break;
	default:
		std::cerr << "Error: unknown sampler" << std::endl;
		return false;
	}

	char line[3000];
	snprintf(line, 3000, "#\n#%14i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i\n", 1,
			2,  3,  4,  5,  6,  7,
			8,  9,  10, 11, 12, 13,
			14, 15, 16, 17, 18, 19,
			20, 21, 22, 23, 24, 25);
	out << line;
	snprintf(line, 3000, "#%14s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n", "time",
			"mean_lq",    "var_lq",    "cnt_lq",    "median_lq",    "min_lq",    "max_lq",
			"mean_gth",   "var_gth",   "cnt_gth",   "median_gth",   "min_gth",   "max_gth",
			"mean_power", "var_power", "cnt_power", "median_power", "min_power", "max_power",
			"mean_other", "var_other", "cnt_other", "median_other", "min_other", "max_other"
			);
	out << line;
	for (size_t i = 0; i < error_acc.size(); i++) {
		size_t time  = error_acc[i].step;
		double mean_lq    = boost::accumulators::mean(error_acc[i].err1);
		double mean_gth   = boost::accumulators::mean(error_acc[i].err2);
		double mean_power = boost::accumulators::mean(error_acc[i].err3);
		double mean_other = boost::accumulators::mean(error_acc[i].err4);
		double var_lq     = boost::accumulators::variance(error_acc[i].err1);
		double var_gth    = boost::accumulators::variance(error_acc[i].err2);
		double var_power  = boost::accumulators::variance(error_acc[i].err3);
		double var_other  = boost::accumulators::variance(error_acc[i].err4);
		size_t cnt_lq     = boost::accumulators::count(error_acc[i].err1);
		size_t cnt_gth    = boost::accumulators::count(error_acc[i].err2);
		size_t cnt_power  = boost::accumulators::count(error_acc[i].err3);
		size_t cnt_other  = boost::accumulators::count(error_acc[i].err4);
		double mdn_lq     = boost::accumulators::median(error_acc[i].err1);
		double mdn_gth    = boost::accumulators::median(error_acc[i].err2);
		double mdn_power  = boost::accumulators::median(error_acc[i].err3);
		double mdn_other  = boost::accumulators::median(error_acc[i].err4);
		double min_lq     = boost::accumulators::min(error_acc[i].err1);
		double min_gth    = boost::accumulators::min(error_acc[i].err2);
		double min_power  = boost::accumulators::min(error_acc[i].err3);
		double min_other  = boost::accumulators::min(error_acc[i].err4);
		double max_lq     = boost::accumulators::max(error_acc[i].err1);
		double max_gth    = boost::accumulators::max(error_acc[i].err2);
		double max_power  = boost::accumulators::max(error_acc[i].err3);
		double max_other  = boost::accumulators::max(error_acc[i].err4);
		//double qm_err   = error_acc[i].get<5>() / runs;
		snprintf(line, 3000, "%15lu%15g%15g%15lu%15g%15g%15g%15g%15g%15lu%15g%15g%15g%15g%15g%15lu%15g%15g%15g%15g%15g%15lu%15g%15g%15g\n", time,
				mean_lq,    var_lq,    cnt_lq,    mdn_lq,    min_lq,    max_lq,
				mean_gth,   var_gth,   cnt_gth,   mdn_gth,   min_gth,   max_gth,
				mean_power, var_power, cnt_power, mdn_power, min_power, max_power,
				mean_other, var_other, cnt_other, mdn_other, min_other, max_other
				);
		out << line;

	}
	return true;
}
