#ifndef STATISTICS_H_0PD5PSM1
#define STATISTICS_H_0PD5PSM1

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "rhab/accumulator.h"

#include <vector>

namespace rhab {

struct StepStatistics {
	size_t step;
	AccumulatorMinMax err[5];
	double wl_f;
};

class ErrorAcc {
private:
	int world_rank;
	int world_size;

	typedef std::vector< StepStatistics > error_vec_t;
	error_vec_t error_acc;

public:
	ErrorAcc() {
#ifdef USE_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
		MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#else
		world_rank = 0;
		world_size = 1;
#endif
	}

	void resize(const size_t &sz) {
		error_acc.resize(sz);
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
		error_acc[idx].err[n](value);
	}

	void pull() {
#ifdef USE_MPI
		if (world_rank == 0) {
			for (int i = 1; i < world_size; i++) {
				MPI_Barrier(MPI_COMM_WORLD);
				for (int j = 0; j < error_acc.size(); j++) {
					for (int k = 0; k < 5; k++) {
						AccumulatorMinMax tmp;
						MPI_Recv((void*)&tmp,
								sizeof(AccumulatorMinMax), MPI_BYTE,
								i, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						error_acc[j].err[k] += tmp;
					}
				}
			}
		} else {
			for (int i = 1; i < world_size; i++) {
				MPI_Barrier(MPI_COMM_WORLD);
				if (i == world_rank) {
					for (int j = 0; j < error_acc.size(); j++) {
						for (int k = 0; k < 5; k++) {
							MPI_Send((void*)&(error_acc[j].err[k]),
									sizeof(AccumulatorMinMax), MPI_BYTE,
									0, 42, MPI_COMM_WORLD);
						}
					}
				}
			}
		}
#endif
	}

	void set_step(const size_t &index, const size_t &step) {
		error_acc[index].step = step;
	}

	void set_parameter(const size_t &index, const double &param) {
		error_acc[index].wl_f = param;
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

