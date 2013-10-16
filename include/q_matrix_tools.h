#ifndef Q_MATRIX_TOOLS_H

#define Q_MATRIX_TOOLS_H

#include "typedefs.h"

#include <boost/tuple/tuple.hpp>

namespace rhab {
	/**
	 * normalize_q
	 */
	void normalize_q(const matrix_int_t & Q, matrix_double_t & Qd);

	/**
	 * normalize
	 */
	void normalize(vector_double_t &vec);

	/**
	 * calculate_dos_gth
	 *
	 * Calculate the Density of States by using the GTH Method.
	 */
	bool calculate_dos_gth(matrix_double_t & inner_mat, vector_double_t &dos);

	/**
	 * calculate_dos_power
	 *
	 * Calculate the Density of States by using the Power Iteration.
	 */
	bool calculate_dos_power(const matrix_double_t &inner_mat, vector_double_t &dos);

	/**
	 * calculate_dos_levmar
	 *
	 * Calculates the Density of States by solving a least squares problem.
	 */
	bool calculate_dos_leastsquares(matrix_int_t imat, matrix_double_t &dmat, vector_double_t &dos);

	/**
	 * calculate_error_q
	 */
	boost::tuple<double, double, double, bool, bool, bool>
	calculate_error_q(const vector_double_t &exact, const matrix_int_t &Q, matrix_double_t &Qd);

	/**
	 * calculate_error
	 */
	double calculate_error(const vector_double_t &exact,
			const vector_double_t &dos,
			bool normalize = false);
}

#endif /* end of include guard: Q_MATRIX_TOOLS_H */
