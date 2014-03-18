#ifndef SIMULATION_SYSTEM_H

#define SIMULATION_SYSTEM_H

// C++ Standard Library
#include <fstream>
#include <numeric>

// Boost Program Options
#include <boost/program_options.hpp>

// Boost Format
#include <boost/format.hpp>

#include "simulation_rng.h"
#include "typedefs.h"
#include "rhab/statistics.h"

class SimulationSystem {
protected:
	//! type of error accumulator array
	typedef rhab::ErrorAcc error_acc_t;

	//! Map of std::string <-> boost::any
	settings_t &settings;

	//! Number of steps to perform per run
	size_t steps;

	//! Number of runs to perform
	size_t runs;

	//! Number of bins used for the different matrices and histograms
	size_t n_bins;

	//! Initial error check frequency
	size_t error_check_f;

	//! Tag to be appended to output file names
	std::string tag;

	//! transition counting matrix
	matrix_int_t Q;

	//! normalized transition matrix
	matrix_double_t Qd;

	//! template for filename
	std::string fn_template;

	//! string stream for header
	std::ostringstream out;

	//! error accumulator array
	error_acc_t error_acc;

	//! error accumulator "matrix" to record per-bin-errors from Least Squares Method
	error_mat_t error_per_bin_lsq;

	//! error accumulator "matrix" to record per-bin-errors from GTH Method
	error_mat_t error_per_bin_gth;

	//! error accumulator "matrix" to record per-bin-errors from Power Method
	error_mat_t error_per_bin_pow;

	//! error accumulator "matrix" to record per-bin-errors from Wang Landau Method
	error_mat_t error_per_bin_wl;

	//! Tuple with pointers to error_per_bin_pow, error_per_bin_gth & error_per_bin_lsq
	error_mat_tuple_t error_matrices;

	//! rank inside MPI_COMM_WORLD
	int world_rank;

	//! size of MPI_COMM_WORLD
	int world_size;

	virtual void setup_output() = 0;
	void open_output_file(std::ofstream& o, const std::string& fn);
	virtual void write_header(std::ofstream& o);

public:
	SimulationSystem(settings_t &s);

	virtual boost::program_options::options_description get_program_options() = 0;
	virtual void parse_arguments(boost::program_options::variables_map &vm) = 0;
	virtual bool run() = 0;
	virtual void setup();
	virtual void write_output(size_t run, const std::string &add);

#ifdef USE_MPI
	void combine_final_dos(matrix_double_t&, matrix_double_t&,
			matrix_double_t&, matrix_double_t&,
			const size_t&, const size_t&);
	void combine_err_matrix(error_mat_tuple_t error_matrices);
#endif

	SimulationRNG rng;

};

#endif /* end of include guard: SIMULATION_SYSTEM_H */
