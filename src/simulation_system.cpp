#include "simulation_system.h"
#include "config.h"
#include "GitSHA1.h"
#include <sstream>
#include <iomanip>

// Boost Format
#include <boost/format.hpp>

SimulationSystem::SimulationSystem(settings_t &s)
	: settings(s), error_matrices(&error_per_bin_lsq, &error_per_bin_gth, &error_per_bin_pow, &error_per_bin_wl) {}


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
	error_per_bin_lsq.resize(error_acc_size, n_bins);
	error_per_bin_gth.resize(error_acc_size, n_bins);
	error_per_bin_pow.resize(error_acc_size, n_bins);

	size_t sampler = boost::any_cast<size_t>(settings["sampler"]);
	if (sampler == 1 || sampler == 5) {
		error_per_bin_wl.resize(error_acc_size, n_bins);
	}

	// safety check
	if (steps % error_check_f != 0) {
		std::cerr << "Error: check-frequency and number of steps don't match" << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

void SimulationSystem::open_output_file(std::ofstream &o, const std::string & fn) {
	o.open(fn.c_str());
	if (!o.good()) {
		throw std::runtime_error("Error: could not open output file '" + fn + "'");
	}
	write_header(o);
}

void SimulationSystem::write_header(std::ofstream& out) {
	out << "# Version: " << VERSION_MAJOR << "." << VERSION_MINOR << std::endl;
	out << "# git SHA: " << g_GIT_SHA1 << std::endl;
	out << "# cmdline: " << boost::any_cast<std::string>(settings["cmdline"]) << std::endl;
	out << "# seed: "    << boost::any_cast<size_t>(settings["seed"]) << std::endl;

}

void write_per_bin_error_file(std::ofstream& out, error_mat_t& error_per_bin, const std::string& header) {
	out << header;
	for (size_t j = 0; j < error_per_bin.size2(); j++) {
		out << std::setw(10) << j;
		for (size_t i = 0; i < error_per_bin.size1(); i++) {
			out << std::setw(12) << error_per_bin(i,j).mean();
		}
		out << "\n";
	}

}

void SimulationSystem::write_output(size_t run, const std::string &add) {
	boost::format fmt(fn_template);
	std::string fn = str( fmt % run );
	//
	//! Output file descriptor for general output
	std::ofstream out_main;

	//! Output file descriptor for per-bin-errors from Least Squares Method
	std::ofstream out_lsq;

	//! Output file descriptor for per-bin-errors from GTH Method
	std::ofstream out_gth;

	//! Output file descriptor for per-bin-errors from Power Method
	std::ofstream out_pow;

	//! Output file descriptor for per-bin-errors from Wang Landau Method
	std::ofstream out_wl;

	open_output_file(out_main, fn);
	open_output_file(out_lsq,  fn + ".lsq");
	open_output_file(out_gth,  fn + ".gth");
	open_output_file(out_pow,  fn + ".pow");

	out_main << out.str();
	out_main << "\n# runs:        " << run;
	out_main << "\n#" << std::endl;

	std::ostringstream oss;
	oss << "# " << std::setw(8) << "bin";
	for (size_t i = 0; i < error_acc.size(); i++) {
		oss << std::setw(12) << (i+1);
	}
	oss << "\n";
	oss << "# " << std::setw(8) << "time->";
	char line[3000];
	snprintf(line, 3000, "#\n#%14i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i%15i\n", 1,
			2,  3,  4,  5,  6,
			7,  8,  9,  10, 11,
			12, 13, 14, 15, 16,
			17, 18, 19, 20, 21,
			22, 23, 24, 25, 26,
			27);
	out_main << line;
	snprintf(line, 3000, "#%14s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n", "time",
			"mean_lq",    "var_lq",    "cnt_lq",    "min_lq",    "max_lq",
			"mean_gth",   "var_gth",   "cnt_gth",   "min_gth",   "max_gth",
			"mean_power", "var_power", "cnt_power", "min_power", "max_power",
			"mean_other", "var_other", "cnt_other", "min_other", "max_other",
			"mean_q",     "var_q",     "cnt_q",     "min_q",     "max_q",
			"wl_f"
			);
	out_main << line;
	for (size_t i = 0; i < error_acc.size(); i++) {
		size_t time  = error_acc[i].step;
		double mean_lq    = boost::accumulators::mean(error_acc[i].err1);
		double mean_gth   = boost::accumulators::mean(error_acc[i].err2);
		double mean_power = boost::accumulators::mean(error_acc[i].err3);
		double mean_other = boost::accumulators::mean(error_acc[i].err4);
		double mean_q     = boost::accumulators::mean(error_acc[i].err_q);
		double var_lq     = boost::accumulators::variance(error_acc[i].err1);
		double var_gth    = boost::accumulators::variance(error_acc[i].err2);
		double var_power  = boost::accumulators::variance(error_acc[i].err3);
		double var_other  = boost::accumulators::variance(error_acc[i].err4);
		double var_q      = boost::accumulators::variance(error_acc[i].err_q);
		size_t cnt_lq     = boost::accumulators::count(error_acc[i].err1);
		size_t cnt_gth    = boost::accumulators::count(error_acc[i].err2);
		size_t cnt_power  = boost::accumulators::count(error_acc[i].err3);
		size_t cnt_other  = boost::accumulators::count(error_acc[i].err4);
		size_t cnt_q      = boost::accumulators::count(error_acc[i].err_q);
		double min_lq     = boost::accumulators::min(error_acc[i].err1);
		double min_gth    = boost::accumulators::min(error_acc[i].err2);
		double min_power  = boost::accumulators::min(error_acc[i].err3);
		double min_other  = boost::accumulators::min(error_acc[i].err4);
		double min_q      = boost::accumulators::min(error_acc[i].err_q);
		double max_lq     = boost::accumulators::max(error_acc[i].err1);
		double max_gth    = boost::accumulators::max(error_acc[i].err2);
		double max_power  = boost::accumulators::max(error_acc[i].err3);
		double max_other  = boost::accumulators::max(error_acc[i].err4);
		double max_q      = boost::accumulators::max(error_acc[i].err_q);
		//double qm_err   = error_acc[i].get<5>() / runs;
		snprintf(line, 3000, "%15lu%15g%15g%15lu%15g%15g%15g%15g%15lu%15g%15g%15g%15g%15lu%15g%15g%15g%15g%15lu%15g%15g%15g%15g%15lu%15g%15g%15g\n", time,
				mean_lq,    var_lq,    cnt_lq,    min_lq,    max_lq,
				mean_gth,   var_gth,   cnt_gth,   min_gth,   max_gth,
				mean_power, var_power, cnt_power, min_power, max_power,
				mean_other, var_other, cnt_other, min_other, max_other,
				mean_q,     var_q,     cnt_q,     min_q,     max_q,
				error_acc[i].wl_f
				);
		out_main << line;
		oss << std::setw(12) << time;
	}

	out_main << add;

	oss << "\n";
	write_per_bin_error_file(out_lsq, error_per_bin_lsq, oss.str());
	write_per_bin_error_file(out_gth, error_per_bin_gth, oss.str());
	write_per_bin_error_file(out_pow, error_per_bin_pow, oss.str());

	size_t sampler = boost::any_cast<size_t>(settings["sampler"]);
	if (sampler == 1 || sampler == 5) {
		write_per_bin_error_file(out_wl,  error_per_bin_wl,  oss.str());
	}

}
