#ifndef SIMULATION_SYSTEM_H

#define SIMULATION_SYSTEM_H

// C++ Standard Library
#include <fstream>

// Boost Program Options
#include <boost/program_options.hpp>

#include "simulation_rng.h"
#include "typedefs.h"
#include "rhab/statistics.h"

class SimulationSystem {
protected:
	// type of error accumulator array
	typedef std::vector< rhab::StepStatistics > error_acc_t;

	settings_t &settings;
	size_t steps;
	size_t runs;

	//! Number of bins used for the different matrices and histograms
	size_t n_bins;

	size_t error_check_f;
	std::ofstream out;
	std::string tag;

	// transition counting matrix
	matrix_int_t Q;

	// normalized transition matrix
	matrix_double_t Qd;

	// error accumulator array
	error_acc_t error_acc;

	//! error accumulator "matrix" to record per-bin-errors from Least Squares Method
	error_mat_t error_per_bin_lsq;

	//! error accumulator "matrix" to record per-bin-errors from GTH Method
	error_mat_t error_per_bin_gth;

	//! error accumulator "matrix" to record per-bin-errors from Power Method
	error_mat_t error_per_bin_pow;

	//! Tuple with pointers to error_per_bin_pow, error_per_bin_gth & error_per_bin_lsq
	error_mat_tuple_t error_matrices;

	virtual void setup_output() = 0;
	bool open_output_files(const std::string& fn);
	bool open_output_file(std::ofstream& o, const std::string& fn);
	virtual void write_header();

public:
	SimulationSystem(settings_t &s);

	virtual boost::program_options::options_description get_program_options() = 0;
	virtual void parse_arguments(boost::program_options::variables_map &vm) = 0;
	virtual bool run() = 0;
	virtual void setup();
	virtual void write_output();

	SimulationRNG rng;

};

#endif /* end of include guard: SIMULATION_SYSTEM_H */
