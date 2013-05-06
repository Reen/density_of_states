#include "q_matrix_tools.h"
#include "rhab/basic_iteration.h"
#include <boost/numeric/ublas/operation.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

matrix_double_t rhab::normalize_q(const matrix_int_t & Q) {
	using namespace boost::numeric::ublas;
	matrix_double_t Qd(Q);
	for (size_t i = 0; i < Qd.size1(); i++) {
		double s = sum(row(Qd,i));
		if (s == 0) continue;
		row(Qd,i) /= s;
	}
	return Qd;
}

void rhab::normalize(vector_double_t &vec) {
	vec /= sum(vec);
}

vector_double_t rhab::calculate_dos_gth(matrix_double_t inner_mat) {
	namespace ublas = boost::numeric::ublas;
	std::size_t inner_rows(inner_mat.size1());
	std::size_t inner_cols(inner_mat.size1());
	vector_double_t dos(inner_mat.size1());
	// we assume small matrix and try GTH method
	// do GTH LU decomposition
	for (std::size_t i = inner_rows-1; i > 0; --i) {
		double s = ublas::norm_1(ublas::subrange(ublas::row(inner_mat,i), 0, i));
		if (s != 0) {
			inner_mat(i,i) = -s;
			for (std::size_t j = 0; j < i; ++j) {
				inner_mat(j,i) /= s;
			}
		}
		for (std::size_t k = 0; k < i; ++k) {
			for (std::size_t j = 0; j < i; ++j) {
				inner_mat(k,j) += inner_mat(k,i)*inner_mat(i,j);
			}
		}
	}
	// now just do eqn 33 of M. Fenwick - J. Chem. Phys. 125, 144905
	std::fill(dos.begin(), dos.end(), 0.0);//dos.clear();
	dos[0] = 1;
	for (std::size_t i = 1; i < inner_rows; ++i) {
		//std::cout << i << " " << dos.data()[i-1] << ": ";
		for (std::size_t j = 0; j < i; ++j) {
			//std::cout << "(" << dos.data()[j] << " " << dos.data()[i-1] << " " << log(whole(j,i)) << ") ";
			dos[i] += exp(dos[j] - dos[i-1] + log(inner_mat(j,i)));
		}
		//std::cout << ": "<< dos.data()[i] << "\n";
		dos[i] = dos[i-1] + log(dos[i]);
	}
	for (std::size_t ei = 0; ei < inner_cols; ++ei) {
		dos(ei) = exp(dos(ei)-5);
	}
	//print_dos("gth", dos, 4);
	normalize(dos);
	return dos;
}

vector_double_t rhab::calculate_dos_power(matrix_double_t mat) {
	namespace ublas = boost::numeric::ublas;
	vector_double_t t1(mat.size1());
	vector_double_t t2(mat.size1());
	vector_double_t t3(mat.size1());
	std::fill(t1.begin(), t1.end(), 1.0/mat.size1());
	size_t max_iter = 10000;
	double lambda, residual, dist;

	rhab::basic_iteration<vector_double_t::iterator> iter(max_iter, 1e-8);
	do {
		++iter;
		noalias(t2) = ublas::prod(t1, mat);
		//ublas::axpy_prod(t1, mat, t2, true);
		t2 /= ublas::norm_1(t2);
		++iter;
		noalias(t1) = ublas::prod(t2, mat);
		//ublas::axpy_prod(t2, mat, t1, true);
		t1 /= ublas::norm_1(t1);
	} while(!iter.converged(t2.begin(), t2.end(), t1.begin(), dist));

	return t1;
}


double rhab::calculate_error(const vector_double_t &exact, const vector_double_t &dos, bool normalize) {
	if (dos.size() == 0 || exact.size() == 0) {
		std::cerr << "exact or calculated density of states vector has zero length!" << std::endl;
		return -1;
	}
	vector_double_t::const_iterator i1 = dos.begin();
	vector_double_t::const_iterator i2 = exact.begin();
	double sum = 0.0;
	if (normalize) {
		double norm = 0;
		double max  = *(std::max_element(dos.begin(), dos.end()));
		for (; i1 != dos.end(); i1++) {
			norm += exp(*i1-max);
#ifdef DEBUG
			if (!boost::math::isfinite(norm) || !boost::math::isfinite(*i1)) {
				std::cerr << __FILE__ << ":" << __LINE__ << " " << norm << " " << *i1 << std::endl;
			}
#endif
		}
		i1 = dos.begin();
		for (; i1 != dos.end(); i1++, i2++) {
			sum += fabs((exp(*i1-max)/norm - *i2) / (*i2));
#ifdef DEBUG
			if (!boost::math::isfinite(sum) || !boost::math::isfinite(*i2)) {
				std::cerr << __FILE__ << ":" << __LINE__ << " " << sum << " " << norm << " " << *i1 <<  " " << *i2 << std::endl;
			}
#endif
		}
		//std::cout << std::setprecision(22) << norm << " " << sum << " " << dos << " " << norm << std::endl;
	} else {
		for (; i1 != dos.end(); i1++, i2++) {
			sum += fabs((*i1 - *i2) / (*i2));
#ifdef DEBUG
			if (!boost::math::isfinite(sum) || !boost::math::isfinite(*i1) || !boost::math::isfinite(*i2)) {
				std::cerr << __FILE__ << ":" << __LINE__ << " " << sum << " " << *i1 << " " << *i2 << std::endl;
			}
#endif
		}
	}
	return sum;
}

double rhab::calculate_error_q(const vector_double_t &exact, const matrix_double_t &Qd) {
	vector_double_t dos(calculate_dos_power(Qd));
	return calculate_error(exact, dos);
}
