#ifndef Q_MATRIX_TOOLS_H

#define Q_MATRIX_TOOLS_H

#include "typedefs.h"


namespace rhab {
	matrix_double_t normalize_q(const matrix_int_t & Q);
	void normalize(vector_double_t &vec);
	vector_double_t calculate_dos_gth(matrix_double_t inner_mat);
	double calculate_error_q(const vector_double_t &exact, const matrix_double_t &Qd);
	double calculate_error(const vector_double_t &exact, const vector_double_t &dos, bool normalize = false);
}

#endif /* end of include guard: Q_MATRIX_TOOLS_H */
