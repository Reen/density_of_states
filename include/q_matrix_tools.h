#ifndef Q_MATRIX_TOOLS_H

#define Q_MATRIX_TOOLS_H

#include "typedefs.h"


namespace rhab {
	/**
	 * normalize_q
	 */
	matrix_double_t normalize_q(const matrix_int_t & Q);

	/**
	 * normalize
	 */
	void normalize(vector_double_t &vec);

	/**
	 * calculate_dos_gth
	 *
	 * Calculate the Density of States by using the GTH Method.
	 */
	vector_double_t calculate_dos_gth(matrix_double_t inner_mat);

	/**
	 * calculate_dos_power
	 *
	 * Calculate the Density of States by using the Power Iteration.
	 */
	vector_double_t calculate_dos_power(const matrix_double_t &inner_mat);

	/**
	 * calculate_error_q
	 */
	double calculate_error_q(const vector_double_t &exact, const matrix_double_t &Qd);

	/**
	 * calculate_error
	 */
	double calculate_error(const vector_double_t &exact,
			const vector_double_t &dos,
			bool normalize = false);
}

#endif /* end of include guard: Q_MATRIX_TOOLS_H */
