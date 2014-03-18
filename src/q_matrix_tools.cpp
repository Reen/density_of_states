#include "q_matrix_tools.h"
#include "rhab/basic_iteration.h"
#include <iomanip>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

// Boost Algorithm
#include <boost/algorithm/minmax_element.hpp>

extern "C"
{
#include <cblas.h>
}

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multimin.h>

typedef boost::tuple<const std::vector<size_t>*,
                     const std::vector<double>*,
                     const std::vector<double>*
                     > params_t;

/**
 * Apply symmetry condition.
 *
 * We can only use pairs of T(i->j) and T(j->i) for the equation system.
 */
void apply_symmetry_conditions(matrix_int_t &imat, std::vector<size_t> &coords) {
  for (size_t i = 0; i < imat.size1(); i++) {
    for (size_t j = 0; j < imat.size2(); ++j) {
      if (i==j) {
        continue;
      }
      if (imat(i,j) > 0 && imat(j,i) > 0) {
        coords.push_back(i);
        coords.push_back(j);
        continue;
      }
      imat(i,j) = imat(j,i) = 0;
    }
  }
}

double variance_function(const gsl_vector * x, void * params) {
  params_t* p = (params_t*)params;
  const std::vector<size_t>* coords  = boost::get<0>(*p);
  const std::vector<double>* sigma   = boost::get<1>(*p);
  const std::vector<double>* precalc = boost::get<2>(*p);

  double err = 0;

  for (size_t k = 0; k < coords->size(); k+=2) {
    const size_t &i = (*coords)[k];
    const size_t &j = (*coords)[k+1];

    double xi(1.0), xj(1.0);
    if (i < x->size) {
      xi = gsl_vector_get(x, i);
    }
    if (j < x->size) {
      xj = gsl_vector_get(x, j);
    }

    double tmp = xi - xj + (*precalc)[k/2];
    err += tmp*tmp / (*sigma)[k/2];
  }

  return err;
}

bool rhab::calculate_dos_minimization(matrix_int_t imat, matrix_double_t &dmat, vector_double_t &dos) {
  // Histogram
  vector_int_t hist(imat.size1());

  // coordinates of entries in the matrix stored one after another i1,j1,i2,j2,...
  std::vector<size_t> coords;
  std::vector<double> sigma;
  std::vector<double> precalc;

  apply_symmetry_conditions(imat, coords);

  // if we have too few entries left in the matrix abort
  if (coords.size()/2 < dos.size()) {
    return false;
  }

  normalize_q(imat, dmat);

  for (size_t i = 0; i < imat.size1(); i++) {
    hist[i] = sum(row(imat,i));
  }

  // we keep the state with the highest energy fixed
  // and thus do minimization only on N-1 states
  size_t num_states = dos.size()-1;

  // Precalculate variables
  sigma.resize(coords.size()/2);
  precalc.resize(coords.size()/2);

  for (size_t k = 0; k < coords.size(); k+=2) {
    const size_t &i = coords[k];
    const size_t &j = coords[k+1];

    precalc[k/2] = log(dmat(i,j) / dmat(j,i));
    sigma[k/2]   = 1./imat(i,j) + 1./hist[i] + 1./imat(j,i) + 1./hist[j];
  }

  // setup the minimization function and fminimizer
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  params_t params(&coords, &sigma, &precalc);

  minex_func.n = num_states;
  minex_func.f = variance_function;
  minex_func.params = &params;

  // Starting point
  x = gsl_vector_alloc(num_states);
  gsl_vector_set_all(x,  1.0);

  // Set initial step sizes to 1
  ss = gsl_vector_alloc(num_states);
  gsl_vector_set_all(ss, 1.0);

  gsl_multimin_fminimizer *s = NULL;
  s = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, num_states);
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);
  int status;
  size_t iter = 0;
  double size;

  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);

    if (status) {
      break;
    }

    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, 1e-6);
    //std::cout << iter << " " << size << " " << status << std::endl;
  } while(status == GSL_CONTINUE /* && iter < 1000000*/);

  std::cout << iter << " " << size;

  for (size_t i = 0; i < s->x->size; i++) {
    dos[i] = gsl_vector_get(s->x, i);
  }
  dos[dos.size()-1] = 1;

  for (size_t i = 0; i < dos.size(); i++) {
    std::cout << " " << dos[i];
  }
  std::cout << std::endl;

  gsl_multimin_fminimizer_free(s);
  gsl_vector_free(x);
  gsl_vector_free(ss);


  return true;
}

bool rhab::calculate_dos_leastsquares(matrix_int_t imat, matrix_double_t &dmat, vector_double_t &dos) {
  vector_int_t hist(imat.size1());

  // coordinates of entries in the matrix stored one after another i1,j1,i2,j2,...
  std::vector<size_t> coords;
  std::vector<size_t> pivot(imat.size1());
  size_t cnt_nnz = 0;

  apply_symmetry_conditions(imat, coords);

  if (coords.size()/2 < dos.size()) {
    return false;
  }

  normalize_q(imat, dmat);

  for (size_t i = 0; i < imat.size1(); i++) {
    hist[cnt_nnz] = sum(row(imat,i));
    if (hist[cnt_nnz] > 0) {
      pivot[i] = cnt_nnz;
      cnt_nnz ++;
    }
  }

  int xn = coords.size()/2;

  gsl_matrix *X, *cov;
  gsl_vector *y, *c;
  double chisq;

  X = gsl_matrix_calloc(xn, cnt_nnz);
  y = gsl_vector_calloc(xn);
  c = gsl_vector_calloc(cnt_nnz);
  cov = gsl_matrix_calloc(cnt_nnz, cnt_nnz);

  for (int k = 0; k < 2*xn; k+=2) {
    const size_t &i = coords[k];
    const size_t &j = coords[k+1];
    gsl_matrix_set(X, k/2, pivot[i],  1.0/sqrt(1./hist[pivot[i]] + 1./hist[pivot[j]] + 1./imat(j,i) + 1./imat(i,j)));
    gsl_matrix_set(X, k/2, pivot[j], -1.0/sqrt(1./hist[pivot[i]] + 1./hist[pivot[j]] + 1./imat(j,i) + 1./imat(i,j)));
    gsl_vector_set(y, k/2, -log( dmat(i,j) / dmat(j,i) )/sqrt(1./hist[pivot[i]] + 1./hist[pivot[j]] + 1./imat(j,i) + 1./imat(i,j)));
  }

  {
    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(xn, cnt_nnz);
    gsl_multifit_linear(X,y,c,cov,&chisq, work);
    gsl_multifit_linear_free(work);
  }

  for (size_t i = 0; i < dos.size(); i++) {
    dos[i] = gsl_vector_get(c, pivot[i]);
  }

  gsl_matrix_free(X);
  gsl_vector_free(y);
  gsl_vector_free(c);
  gsl_matrix_free(cov);

  return true;
}

void rhab::normalize_q(const matrix_int_t & Q, matrix_double_t & Qd) {
  using namespace boost::numeric::ublas;
  //matrix_double_t Qd(Q);
  for (size_t i = 0; i < Q.size1(); i++) {
    double s = sum(row(Q,i));
    if (s == 0) continue;
    row(Qd,i) = row(Q,i)/s;
  }
}

void rhab::normalize(vector_double_t &vec) {
  vec /= sum(vec);
}

void rhab::normalize_from_log(vector_double_t &vec) {
  double sub(0), norm(0);
  std::pair<vector_double_t::iterator, vector_double_t::iterator> mm =
    boost::minmax_element(vec.begin(), vec.end());
  sub = (*(mm.second) + *(mm.first))/2;
  for (size_t i = 0; i < vec.size(); i++) {
    norm += exp(vec[i] - sub);
  }
  for (size_t i = 0; i < vec.size(); i++) {
    vec[i] = exp(vec[i] - sub) / norm;
  }
}

bool rhab::calculate_dos_gth(matrix_double_t & inner_mat, vector_double_t &dos) {
  namespace ublas = boost::numeric::ublas;
  std::size_t inner_rows(inner_mat.size1());
  std::size_t inner_cols(inner_mat.size1());
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
  // now just do modified eqn. 33 of M. Fenwick - J. Chem. Phys. 125, 144905
  std::fill(dos.begin(), dos.end(), 0.0);//dos.clear();
  dos[0] = 1;
  for (std::size_t i = 1; i < inner_rows; ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      if (inner_mat(j,i) > 0) {
        dos[i] += exp(dos[j] + log(inner_mat(j,i)));
      }
    }
    dos[i] = log(dos[i]);
  }
  for (std::size_t ei = 0; ei < inner_cols; ++ei) {
    dos(ei) = exp(dos(ei));
  }

  // eqn. 32 of M. Fenwick - J. Chem. Phys. 125, 144905
  /*
  for (size_t i = 1; i < inner_rows; ++i) {
    for (size_t j = 0; j < i; ++j) {
      dos[i] += dos[j] * inner_mat(j,i);
    }
  }
  */

  normalize(dos);
  return ( std::count_if(dos.begin(), dos.end(), boost::math::isnan<double>) == 0 );
}

bool rhab::calculate_dos_power(const matrix_double_t &imat, vector_double_t &t1) {
  namespace ublas = boost::numeric::ublas;
  matrix_double_t mat(ublas::trans(imat));
  //vector_double_t t1(mat.size1());
  vector_double_t t2(mat.size1());
  vector_double_t t3(mat.size1());
  std::fill(t1.begin(), t1.end(), 1.0/mat.size1());
  size_t max_iter = 10000;
  double lambda, residual, dist;

  rhab::basic_iteration<vector_double_t::iterator> iter(max_iter, 1e-8);
  do {
    ++iter;
    cblas_dgemv(CblasRowMajor, CblasNoTrans, mat.size1(), mat.size2(), 1, &(mat.data()[0]), mat.size1(), &(t1.data()[0]), 1, 0.0, &(t2.data()[0]), 1);
    //ublas::axpy_prod(mat, t1, t2, true);
    t2 /= ublas::norm_1(t2);
    ++iter;
    cblas_dgemv(CblasRowMajor, CblasNoTrans, mat.size1(), mat.size2(), 1, &(mat.data()[0]), mat.size1(), &(t2.data()[0]), 1, 0.0, &(t1.data()[0]), 1);
    //ublas::axpy_prod(mat, t2, t1, true);
    t1 /= ublas::norm_1(t1);
  } while(!iter.converged(t2.begin(), t2.end(), t1.begin(), dist));

  size_t cnt_zero = std::count_if(t2.begin(), t2.end(), std::bind2nd(std::equal_to<double>(), 0.0));
  bool enought_non_zero   = ( cnt_zero < 0.2 * t2.size() );
  bool all_entries_finite = ( std::count_if(t2.begin(), t2.end(),
        boost::math::isfinite<double>) == (long)t2.size() );
  return ( enought_non_zero && all_entries_finite );
}


double rhab::calculate_error(const vector_double_t &exact,
                             const vector_double_t &dos,
                             error_mat_t* error_per_bin,
                             const size_t& index, bool normalize) {
  if (dos.size() == 0 || exact.size() == 0) {
    std::cerr << "exact or calculated density of states vector has zero length!" << std::endl;
    return -1;
  }

  // vector to calculate the error
  vector_double_t err(dos.size());
  std::fill(err.begin(), err.end(), 0.0);

  double sum = 0.0;
  size_t cnt = 0;
  double norm = 0;

  if (normalize) {
    /*
     * Density of states provided in dos is actually \f$\ln(\Omega)\f$.
     * Find the largest value, then subtract it from dos[i]
     * and sum exp(dos[i] - max) to calculate the norm.
     * Then divide the every exp(dos[i] - max) by the norm
     * and subtract the exact value, i.e. exact[i].
     * Calculate the absolute value of it and divide by exact[i].
     *
     * exact[i] is assumed to be positive
     *
     * @todo: calculate both error in density of states and entropy,
     *        i.e. additionally dos[i]-log(exact[i])/log(exact[i])
     */

    vector_double_t d(dos);
    vector_double_t::iterator middle = d.begin()+(d.end()-d.begin())/2;
    std::nth_element(d.begin(), middle, d.end());
    double min(0), max(0);
    double sub  = *middle;

    // calculate the norm
    for (size_t i = 0; i < dos.size(); i++) {
      // Be careful here! The least squares and minimization algorithms
      // tend to calculate very large values for states that have never been
      // visited. Having a very large max value results in the exp(foo-max)
      // to become 0
      if ((exact[i]) > 0) {
        norm += exp(dos[i]-sub);
#ifdef DEBUG
        if (!boost::math::isfinite(norm) || !boost::math::isfinite(dos[i])
            || !boost::math::isfinite(exp(dos[i]-sub))) {
          std::cerr << __FILE__ << ":" << __LINE__ << " "
                    << norm << " " << dos[i] << " " << exp(dos[i]-sub) << std::endl;
        }
#endif
      }
    }

    // yeah, I know that one should not compare doubles by equal
    if (norm == 0) {
      norm = 1;
    }

    for (size_t i = 0; i < dos.size(); i++) {
      // Be careful here and do not devide by 0
      if ((exact[i]) > 0) {
        err[i] = fabs( (exp(dos[i]-sub)/norm - exact[i]) / exact[i] );
        sum += err[i];
        cnt ++;
        (*error_per_bin)(index, i)(err[i]);
#ifdef DEBUG
        if (!boost::math::isfinite(sum) || !boost::math::isfinite(exact[i])
            || !boost::math::isfinite(err[i])) {
          std::cerr << __FILE__ << ":" << __LINE__ << " "
                    << sum << " " << norm << " "
                    << dos[i] <<  " " << exact[i] << std::endl;
        }
#endif
      }
    }
  } else {
    for (size_t i = 0; i < dos.size(); i++) {
      if (exact[i] > 0) {
        norm += dos[i];
      }
    }
    for (size_t i = 0; i < dos.size(); i++) {
      // Be careful here and do not divide by 0
      if (exact[i] > 0) {
        err[i] = fabs( (dos[i]/norm - exact[i]) / exact[i] );
        sum += err[i];
        cnt ++;
        (*error_per_bin)(index, i)(err[i]);
      }
#ifdef DEBUG
      if (!boost::math::isfinite(sum) || !boost::math::isfinite(dos[i])
          || !boost::math::isfinite(exact[i])
          || !boost::math::isfinite(err[i])) {
        std::cerr << __FILE__ << ":" << __LINE__ << " "
                  << sum << " " << dos[i] << " " << exact[i] << std::endl;
      }
#endif
    }
  }

  return (sum / cnt);
}

double rhab::calculate_error_q_matrix(const matrix_double_t &Qex, const matrix_double_t &Q) {
  double value(0.0);
  size_t cnt(0);
  for (size_t i = 0; i < Qex.data().size(); i++) {
    if (Qex.data()[i] > 0) {
      double tmp = fabs(Qex.data()[i] - Q.data()[i])/Qex.data()[i];
      value += tmp;
      cnt ++;
    }
  }
  return value / cnt;
}

boost::tuple<double, double, double, bool, bool, bool, double>
rhab::calculate_error_q(const vector_double_t &exact, const matrix_double_t &Qexact,
                        const matrix_int_t &Q, matrix_double_t &Qd,
                        error_mat_tuple_t error_matrices,
                        const size_t& index) {
  vector_double_t dos(Q.size1());
  std::fill(dos.begin(), dos.end(), 0.0);

  // Least Squares
  bool lq = calculate_dos_leastsquares(Q, Qd, dos);
  //bool lq = calculate_dos_minimization(Q, Qd, dos);
  double error_lsq = calculate_error(exact, dos, error_matrices.get<0>(), index, true);
  if (!boost::math::isfinite(error_lsq)) {
    lq = false;
  }

  // Least Squares uses Qd as workspace only, so compute Qd
  normalize_q(Q, Qd);
  std::fill(dos.begin(), dos.end(), 0.0);

  // calculate error in Q matrix
  double q_error = 0;
  if (Qexact.size1() > 0) {
    q_error = rhab::calculate_error_q_matrix(Qexact, Qd);
  }

  // GTH method
  bool gth = calculate_dos_gth(Qd, dos);
  double error_gth = calculate_error(exact, dos, error_matrices.get<1>(), index, false);

  // GTH Method modifies Qd, so recompute
  normalize_q(Q, Qd);
  std::fill(dos.begin(), dos.end(), 0.0);

  // Power method
  bool pow = calculate_dos_power(Qd, dos);
  double error_pow = calculate_error(exact, dos, error_matrices.get<2>(), index, false);

  return boost::make_tuple(error_lsq, error_gth, error_pow, lq, gth, pow, q_error);
}

boost::tuple<double, double, double, bool,  bool,   bool>
rhab::calculate_error_q_lj(const vector_double_t &exact,
                           const matrix_int_t &Q, matrix_double_t &Qd,
                           error_mat_tuple_t error_matrices,
                           vector_double_t &dos_lsq,
                           vector_double_t &dos_gth,
                           vector_double_t &dos_pow,
                           const size_t& index) {
  std::fill(dos_lsq.begin(), dos_lsq.end(), 0.0);
  std::fill(dos_gth.begin(), dos_gth.end(), 0.0);
  std::fill(dos_pow.begin(), dos_pow.end(), 0.0);

  // Least Squares
  bool lq = calculate_dos_leastsquares(Q, Qd, dos_lsq);
  //bool lq = calculate_dos_minimization(Q, Qd, dos_lsq);
  double error_lsq = calculate_error(exact, dos_lsq, error_matrices.get<0>(), index, true);
  if (!boost::math::isfinite(error_lsq)) {
    lq = false;
  }

  // Least Squares uses Qd as workspace only, so compute Qd
  normalize_q(Q, Qd);

  // GTH method
  bool gth = calculate_dos_gth(Qd, dos_gth);
  double error_gth = calculate_error(exact, dos_gth, error_matrices.get<1>(), index, false);

  // GTH Method modifies Qd, so recompute
  normalize_q(Q, Qd);

  // Power method
  bool pow = calculate_dos_power(Qd, dos_pow);
  double error_pow = calculate_error(exact, dos_pow, error_matrices.get<2>(), index, false);

  return boost::make_tuple(error_lsq, error_gth, error_pow,
                           lq,        gth,       pow);
}


/* vim: set ts=2 sw=2 sts=2 tw=0 expandtab :*/
