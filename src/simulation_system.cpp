#include "simulation_system.h"
#include "config.h"
#include "GitSHA1.h"

SimulationSystem::SimulationSystem(settings_t &s) : settings(s) {}


void SimulationSystem::setup() {
	if (settings.count("seed")) {
		rng.seed = boost::any_cast<size_t>(settings["seed"]);
		rng.seed_set = true;
	} else {
		size_t seed = rng.initialize();
		settings["seed"] = seed;
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

	size_t error_acc_size = 9 * (size_t)log10(steps / error_check_f) + 1;
	error_acc.resize(error_acc_size);

	// safety check
	if (steps % error_check_f != 0) {
		std::cerr << "Error: check-frequency and number of steps don't match" << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

bool SimulationSystem::open_output_files(const std::string& fn) {
	open_output_file(out,     fn);
}

bool SimulationSystem::open_output_file(std::ofstream &o, const std::string & fn) {
	o.open(fn.c_str());
	if (!o.good()) {
		throw std::runtime_error("Error: could not open output file '" + fn + "'");
	}
	return true;
}

void SimulationSystem::write_header() {
	out << "# Version: " << VERSION_MAJOR << "." << VERSION_MINOR << std::endl;
	out << "# git SHA: " << g_GIT_SHA1 << std::endl;
	out << "# cmdline: " << boost::any_cast<std::string>(settings["cmdline"]) << std::endl;
	out << "# seed: "    << boost::any_cast<size_t>(settings["seed"]) << std::endl;
}

void SimulationSystem::write_output() {
	char line[3000];
	snprintf(line, 3000, "#\n#%14i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i\n", 1,
			2,  3,  4,  5,  6,
			7,  8,  9,  10, 11,
			12, 13, 14, 15, 16,
			17, 18, 19, 20, 21);
	out << line;
	snprintf(line, 3000, "#%14s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n", "time",
			"mean_lq",    "var_lq",    "cnt_lq",    "min_lq",    "max_lq",
			"mean_gth",   "var_gth",   "cnt_gth",   "min_gth",   "max_gth",
			"mean_power", "var_power", "cnt_power", "min_power", "max_power",
			"mean_other", "var_other", "cnt_other", "min_other", "max_other"
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
		double min_lq     = boost::accumulators::min(error_acc[i].err1);
		double min_gth    = boost::accumulators::min(error_acc[i].err2);
		double min_power  = boost::accumulators::min(error_acc[i].err3);
		double min_other  = boost::accumulators::min(error_acc[i].err4);
		double max_lq     = boost::accumulators::max(error_acc[i].err1);
		double max_gth    = boost::accumulators::max(error_acc[i].err2);
		double max_power  = boost::accumulators::max(error_acc[i].err3);
		double max_other  = boost::accumulators::max(error_acc[i].err4);
		//double qm_err   = error_acc[i].get<5>() / runs;
		snprintf(line, 3000, "%15lu%15g%15g%15lu%15g%15g%15g%15g%15lu%15g%15g%15g%15g%15lu%15g%15g%15g%15g%15lu%15g%15g\n", time,
				mean_lq,    var_lq,    cnt_lq,    min_lq,    max_lq,
				mean_gth,   var_gth,   cnt_gth,   min_gth,   max_gth,
				mean_power, var_power, cnt_power, min_power, max_power,
				mean_other, var_other, cnt_other, min_other, max_other
				);
		out << line;

	}
}
