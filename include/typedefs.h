#ifndef TYPEDEFS_H

#define TYPEDEFS_H

// C++ Standard Library
#include <map>
#include <string>

// Boost Numeric Ublas
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

typedef boost::numeric::ublas::vector<int> vector_int_t;
typedef boost::numeric::ublas::vector<double> vector_double_t;
typedef boost::numeric::ublas::matrix<int64_t> matrix_int_t;
typedef boost::numeric::ublas::matrix<double> matrix_double_t;

// Boost Any
#include <boost/any.hpp>

typedef std::map< std::string, boost::any > settings_t;

#include "rhab/accumulator.h"
#include <boost/tuple/tuple.hpp>
typedef boost::numeric::ublas::matrix<rhab::Accumulator> error_mat_t;
typedef boost::tuple<
	error_mat_t*,
	error_mat_t*,
	error_mat_t*,
	error_mat_t*
> error_mat_tuple_t;

#endif /* end of include guard: TYPEDEFS_H */
