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
 * normalize_from_log
 */
void normalize_from_log(vector_double_t &vec);

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
 * calculate_dos_minimization
 *
 * Calculates the Density of States by minimizing a function
 */
bool calculate_dos_minimization(matrix_int_t imat, matrix_double_t &dmat, vector_double_t &dos);

/**
 * calculate_error_q_matrix
 */
double calculate_error_q_matrix(const matrix_double_t &Qex, const matrix_double_t &Q);

/**
 * calculate_error_q
 */
boost::tuple<double, double, double, bool, bool, bool, double>
calculate_error_q(const vector_double_t &exact,
                  const matrix_double_t &Qexact,
                  const matrix_int_t &Q, matrix_double_t &Qd,
                  error_mat_tuple_t error_matrices,
                  vector_double_t &dos_lsq,
                  vector_double_t &dos_gth,
                  vector_double_t &dos_pow,
                  const size_t& index);

/**
 * calculate_error_q_lj
 *
 * Error calculation for the LJ system.
 */
boost::tuple<double, double, double, bool,  bool,   bool>
calculate_error_q_lj(const vector_double_t &exact,
                     const matrix_int_t &Q, matrix_double_t &Qd,
                     error_mat_tuple_t error_matrices,
                     vector_double_t &dos_lsq,
                     vector_double_t &dos_gth,
                     vector_double_t &dos_pow,
                     const size_t& index);

/**
 * calculate_error
 */
double calculate_error(const vector_double_t &exact,
                       const vector_double_t &dos,
                       error_mat_t* error_per_bin, const size_t& index,
                       bool normalize = false);
}

#endif /* end of include guard: Q_MATRIX_TOOLS_H */

/* vim: set ts=2 sw=2 sts=2 tw=0 expandtab :*/
