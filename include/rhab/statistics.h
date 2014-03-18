#ifndef STATISTICS_H_0PD5PSM1
#define STATISTICS_H_0PD5PSM1

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <vector>

// Boost Accumulator
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

namespace rhab {

struct StepStatistics {
	size_t step;
	boost::accumulators::accumulator_set<
		double,
		boost::accumulators::stats<
			boost::accumulators::tag::min,
			boost::accumulators::tag::max,
			boost::accumulators::tag::mean,
			boost::accumulators::tag::variance
				> > err[5];
	double wl_f;
};

class ErrorAcc {
private:
	int world_rank;
	int world_size;

	typedef std::vector< StepStatistics > error_vec_t;
	error_vec_t error_acc;
#ifdef USE_MPI
	MPI_Request rqst;
	MPI_Status stts;
#endif

public:
	ErrorAcc() {
#ifdef USE_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
		MPI_Comm_size(MPI_COMM_WORLD, &world_size);
		rqst = MPI_REQUEST_NULL;
#else
		world_rank = 0;
		world_size = 1;
#endif
	}

	void resize(const size_t &sz) {
		if (world_rank == 0) {
			error_acc.resize(sz);
		}
	}

	struct ErrEntry {
		size_t n;
		size_t idx;
		double value;

		ErrEntry(){}

		ErrEntry(const size_t &n_, const size_t &idx_, const double &value_)
			: n(n_), idx(idx_), value(value_) {}
	};

	// syntax err_acc.push(1, index, value)
	void push(const size_t &n, const size_t &idx, const double &value) {
		if (world_rank == 0) {
			error_acc[idx].err[n](value);
			pull();
		} else {
#ifdef USE_MPI
			ErrEntry ee(n, idx, value);
			// wait for previous request to finish
			//MPI_Wait(&rqst, &stts);
			// perform new non-blocking send
			MPI_Wait(&rqst, &stts);
			MPI_Isend(&ee, sizeof(ErrEntry), MPI_BYTE, 0, 42, MPI_COMM_WORLD, &rqst);
#endif
		}
	}

	void pull() {
#ifdef USE_MPI
		if (world_rank == 0) {
			int flag;
			do {
				MPI_Iprobe(MPI_ANY_SOURCE, 42, MPI_COMM_WORLD, &flag, &stts);
				if (flag) {
					ErrEntry ee;
					MPI_Request rrqst;
					MPI_Irecv(&ee, sizeof(ErrEntry), MPI_BYTE, stts.MPI_SOURCE, 42, MPI_COMM_WORLD, &rrqst);
					MPI_Wait(&rrqst, MPI_STATUS_IGNORE);
					error_acc[ee.idx].err[ee.n](ee.value);
				} else {
					break;
				}
			} while(true);
		}
#endif
	}

	void set_step(const size_t &index, const size_t &step) {
		if (world_rank == 0) {
			error_acc[index].step = step;
		}
	}

	void set_parameter(const size_t &index, const double &param) {
		if (world_rank == 0) {
			error_acc[index].wl_f = param;
		}
	}

	error_vec_t::size_type size() const {
		return error_acc.size();
	}

	StepStatistics& operator[](const size_t& i) {
		return error_acc[i];
	}

};

}

#endif /* end of include guard: STATISTICS_H_0PD5PSM1 */

