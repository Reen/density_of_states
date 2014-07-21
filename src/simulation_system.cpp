#include "simulation_system.h"
#include "config.h"
#include "GitSHA1.h"
#include <sstream>
#include <iomanip>

#ifdef USE_MPI
#include <mpi.h>
#endif

SimulationSystem::SimulationSystem(settings_t &s)
	: settings(s), error_matrices(&error_per_bin_lsq, &error_per_bin_gth, &error_per_bin_pow, &error_per_bin_wl) {
#ifdef USE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#else
	world_rank = 0;
	world_size = 1;
#endif
}


void SimulationSystem::setup() {
	if (settings.count("seed")) {
		rng.seed = boost::any_cast<size_t>(settings["seed"]);
		rng.seed_set = true;
	} else {
		size_t seed = rng.initialize();
		settings["seed"] = seed;
#ifdef USE_MPI
		std::vector<size_t> seedlist(world_size, 0);
		MPI_Barrier(MPI_COMM_WORLD);
		if (world_rank == 0) {
			seedlist[0] = rng.seed;
			for (size_t i = 1; i < (size_t)world_size; i++) {
				MPI_Recv(&seedlist[i], sizeof(size_t), MPI_BYTE, i, 12345, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			for (size_t i = 0; i < (size_t)world_size; i++) {
				seed_out << "# seed[" << i << "] " << seedlist[i] << "\n";
			}
		} else {
			MPI_Send(&rng.seed, sizeof(size_t), MPI_BYTE, 0, 12345, MPI_COMM_WORLD);
		}
#endif
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
	out << seed_out.str();

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
		double mean_lq    = error_acc[i].err[0].mean();
		double mean_gth   = error_acc[i].err[1].mean();
		double mean_power = error_acc[i].err[2].mean();
		double mean_other = error_acc[i].err[3].mean();
		double mean_q     = error_acc[i].err[4].mean();
		double var_lq     = error_acc[i].err[0].variance();
		double var_gth    = error_acc[i].err[1].variance();
		double var_power  = error_acc[i].err[2].variance();
		double var_other  = error_acc[i].err[3].variance();
		double var_q      = error_acc[i].err[4].variance();
		size_t cnt_lq     = error_acc[i].err[0].count();
		size_t cnt_gth    = error_acc[i].err[1].count();
		size_t cnt_power  = error_acc[i].err[2].count();
		size_t cnt_other  = error_acc[i].err[3].count();
		size_t cnt_q      = error_acc[i].err[4].count();
		double min_lq     = error_acc[i].err[0].min();
		double min_gth    = error_acc[i].err[1].min();
		double min_power  = error_acc[i].err[2].min();
		double min_other  = error_acc[i].err[3].min();
		double min_q      = error_acc[i].err[4].min();
		double max_lq     = error_acc[i].err[0].max();
		double max_gth    = error_acc[i].err[1].max();
		double max_power  = error_acc[i].err[2].max();
		double max_other  = error_acc[i].err[3].max();
		double max_q      = error_acc[i].err[4].max();
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

	open_output_file(out_lsq,  fn + ".lsq");
	open_output_file(out_gth,  fn + ".gth");
	open_output_file(out_pow,  fn + ".pow");
	write_per_bin_error_file(out_lsq, error_per_bin_lsq, oss.str());
	write_per_bin_error_file(out_gth, error_per_bin_gth, oss.str());
	write_per_bin_error_file(out_pow, error_per_bin_pow, oss.str());

	size_t sampler = boost::any_cast<size_t>(settings["sampler"]);
	if (sampler == 1 || sampler == 5) {
		open_output_file(out_wl,  fn + ".wl");
		write_per_bin_error_file(out_wl,  error_per_bin_wl,  oss.str());
	}

}



#ifdef USE_MPI

void transfer_matrix(matrix_double_t& fd, const size_t& run_from, const size_t& run_to) {
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (rank == 0) {
		for (int i = 1; i < size; i++) {
			vector_double_t recvbuf(fd.size2());
			MPI_Status stts;
			for (size_t j = run_from; j < run_to; j++) {
				// j is not used here, we just need do do the following
				// for (run_to-run_from) times
				MPI_Recv((void*)recvbuf.data().begin(), recvbuf.size(),
						MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, &stts);
				// the tag transports the acutal run
				//std::cerr << "[" << rank << "] " << i << " " << j << " " << stts.MPI_TAG << " " << stts.MPI_SOURCE << std::endl;
				boost::numeric::ublas::row(fd, stts.MPI_TAG) = recvbuf;
			}
		}
	} else {
		for (size_t r = run_from; r < run_to; r++) {
			//std::cerr << "[" << rank << "] " << r << std::endl;
			vector_double_t sendbuf(boost::numeric::ublas::row(fd, r));
			MPI_Send((void*)sendbuf.data().begin(), sendbuf.size(),
					MPI_DOUBLE, 0, r, MPI_COMM_WORLD);
		}

	}
}

void SimulationSystem::combine_final_dos(matrix_double_t& fd_lsq,
		matrix_double_t& fd_gth, matrix_double_t& fd_pow,
		matrix_double_t& fd_wl,
		const size_t& run_from, const size_t& run_to) {
	//std::cerr << "[" << world_rank << "] in:" << __FUNCTION__ <<":"<<__LINE__ << "\n";
	transfer_matrix(fd_lsq, run_from, run_to);
	MPI_Barrier(MPI_COMM_WORLD);
	transfer_matrix(fd_gth, run_from, run_to);
	transfer_matrix(fd_pow, run_from, run_to);
	transfer_matrix(fd_wl,  run_from, run_to);
	//std::cerr << "[" << world_rank << "] out:" << __FUNCTION__ <<":"<<__LINE__ << std::endl;
}

void transfer_err_matrix(error_mat_t* mat) {
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Barrier(MPI_COMM_WORLD);

	if (mat->size1() == 0 || mat->size2() == 0) {
		return;
	}

	if (rank == 0) {
		error_mat_t tmp(mat->size1(), mat->size2());
		MPI_Status stts;
		std::cout << "err mat size: " << tmp.size1() << " " << tmp.size2() << " " << sizeof(rhab::Accumulator) << std::endl;
		for (int k = 1; k < size; k++) {
			for (size_t i = 0; i < tmp.size1(); i++) {
				for (size_t j = 0; j < tmp.size2(); j++) {
					MPI_Recv((void*)(&tmp(i,j)), sizeof(rhab::Accumulator),
							MPI_BYTE, k, 13, MPI_COMM_WORLD, &stts);
				}
			}
			// add the error matrix of every other node to rank 0's
			(*mat) += tmp;
		}
	} else {
		for (size_t i = 0; i < mat->size1(); i++) {
			for (size_t j = 0; j < mat->size2(); j++) {
				MPI_Send((void*)(&(*mat)(i,j)), sizeof(rhab::Accumulator),
						MPI_BYTE, 0, 13, MPI_COMM_WORLD);
			}
		}
		std::cout << "ERR mat size: " << mat->size1() << " " << mat->size2() << " " << sizeof(rhab::Accumulator) << std::endl;
	}
}

void SimulationSystem::combine_err_matrix(error_mat_tuple_t error_matrices) {
	MPI_Barrier(MPI_COMM_WORLD);
	transfer_err_matrix(error_matrices.get<0>());
	transfer_err_matrix(error_matrices.get<1>());
	transfer_err_matrix(error_matrices.get<2>());
	transfer_err_matrix(error_matrices.get<3>());
}
#endif

