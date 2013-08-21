#ifndef MC_SAMPLER_H

#define MC_SAMPLER_H


// Boost Random
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_01.hpp>

#include "typedefs.h"

/**
 * MCSampler
 */
class MCSampler {
protected:
	boost::mt19937 &rng;
	boost::uniform_01<> dist01;

	settings_t settings;
public:
	MCSampler(boost::mt19937 &_rng, const settings_t &_settings)
		: rng(_rng), settings(_settings) {}

	void check(const size_t & step, const size_t &run) {}
	bool has_own_statistics() {
		return false;
	}
	double calculate_error(const vector_double_t &exact) {
		return 0.0;
	}

	template<class T2>
	void rejected(const T2 &) {}
};


/**
 * BoltzmannFunctor
 */
class BoltzmannFunctor {
	const double &kB;
	const double &T;

public:
	BoltzmannFunctor(const double &_kB, const double &_T)
		: kB(_kB), T(_T) {}

	template<class FloatType>
	inline double operator()(const FloatType &dE) {
		return exp(dE);
	}
};

/**
 * BoltzmannTableFunctor
 */
class BoltzmannTableFunctor {
	const double &kB;
	const double &T;
	double table[9];

public:
	BoltzmannTableFunctor(const double &_kB, const double &_T)
		: kB(_kB), T(_T) {
		table[4] = exp(-4/(kB*T));
		table[8] = exp(-8/(kB*T));
	}

	inline double operator()(const int &dE) {
		BOOST_ASSERT((dE == 4) || (dE == 8));
		return table[dE];
	}
};

/**
 * BoltzmannSampler
 */
template<class expfn_t>
class BoltzmannSampler : public MCSampler {
private:
	double kB;
	double T;
	expfn_t exp_fn;
public:
	BoltzmannSampler(boost::mt19937 &rng, const matrix_int_t&, const settings_t &settings)
		: MCSampler(rng, settings), exp_fn(kB, T) {
		kB = boost::any_cast<double>(settings.find("kB")->second);
		T  = boost::any_cast<double>(settings.find("temperature")->second);
	}

	template<class T1>
	bool operator()(const T1 &dE) {
		return ((dE <= 0) || (dist01(rng) <= exp_fn(dE)));
	}

	template<class T1, class T2>
	bool operator()(const T1 &dE, const T2 &i_old, const T2 &i_new) {
		return ((dE <= 0) || (dist01(rng) <= exp_fn(dE)));
	}

	void check(const size_t & step, const size_t &run) {
		//T = 2 - (1.9 * step/1e7);
		//T = sin(step/2e5)+1.1;
	}
};

/**
 * WangLandauSampler
 */
class WangLandauSampler : public MCSampler {
private:
	vector_int_t H;
	vector_double_t g;
	double ln_f;
	double flatness;
	size_t macrostates;
public:
	WangLandauSampler(boost::mt19937 &rng, const matrix_int_t&, const settings_t &settings)
		: MCSampler(rng, settings), ln_f(1.0)
	{
		macrostates = boost::any_cast<size_t>(settings.find("macrostates")->second);
		flatness    = boost::any_cast<double>(settings.find("flatness")->second);

		H.resize(macrostates);
		g.resize(macrostates);

		H *= 0;
		g *= 0;

		for (size_t i = 0; i < H.size(); i++) {
			assert(H[i] == 0);
			assert(g[i] == 0);
		}
	}

	template<class T1, class T2>
	bool operator()(const T1 &dE, const T2 &i_old, const T2 &i_new) {
		bool res = ((g[i_old] >= g[i_new]) || (dist01(rng) <= exp(g[i_old]-g[i_new])));

		if (res) {
			H[i_new]++;
			g[i_new]+=ln_f;
		} else {
			H[i_old]++;
			g[i_old]+=ln_f;
		}
		return res;
	}

	template<class T2>
	void rejected(const T2 &i_old) {
		H[i_old]++;
		g[i_old]+=ln_f;
	}

	void check(const size_t & step, const size_t &run) {
		int min = *std::min_element(H.begin(), H.end());
		if (min > 0 && (min > (0.99 * sum(H))/macrostates)) {
			H *= 0;
			ln_f /= 2.0;
#if VERBOSE == 1
			std::cout << std::setw(15) << std::right << run
					  << std::setw(15) << std::right << step
					  << std::setw(15) << std::right << ln_f << "\n";
#endif
		}
	}

	bool has_own_statistics() {
		return true;
		//return false;
	}

	double calculate_error(const vector_double_t &exact) {
		return ::calculate_error(exact, g, true);
	}
};

/**
 * QualityMeasureASampler
 */
class QualityMeasureASampler : public MCSampler {
private:
	const matrix_int_t& Q;
public:
	QualityMeasureASampler(boost::mt19937 &rng, const matrix_int_t& qmat, const settings_t &settings)
		: MCSampler(rng, settings), Q(qmat) {
	}

	template<class T1, class T2>
	bool operator()(const T1 &dE, const T2 &i_old, const T2 &i_new) {
		size_t Hold_sum(0), Hnew_sum(0);
		for (size_t i = 0; i < Q.size1(); i++) {
			Hold_sum += Q(i_old,i);
			Hnew_sum += Q(i_new,i);
		}
		bool res = ((Hold_sum >= Hnew_sum) || (dist01(rng) <= exp(Hold_sum-Hnew_sum)));
#if VERBOSE == 1
		//std::cout << i_old << " " << i_new << " " << Hold << " " << Hnew << " " << res << std::endl;
#endif
		return res;
	}

	void check(const size_t &step, const size_t &run) {}

	bool has_own_statistics() {
		return false;
	}
};

/**
 * QualityMeasureBSampler
 */
class QualityMeasureBSampler : public MCSampler {
private:
	const matrix_int_t& Q;
public:
	QualityMeasureBSampler(boost::mt19937 &rng, const matrix_int_t& qmat, const settings_t &settings)
		: MCSampler(rng, settings), Q(qmat) {
	}

	template<class T1, class T2>
	bool operator()(const T1 &dE, const T2 &i_old, const T2 &i_new) {
		size_t Hold_sum(0), Hnew_sum(0), Hold_cnt(0), Hnew_cnt(0);
		for (size_t i = 0; i < Q.size1(); i++) {
			if (Q(i_old,i) != 0) {
				Hold_sum += Q(i_old,i);
				Hold_cnt += 1;
			}
			if (Q(i_new,i) != 0) {
				Hnew_sum += Q(i_new,i);
				Hnew_cnt += 1;
			}
		}
		double Hold = 0;
		double Hnew = 0;
		if (Hold_cnt != 0) {
			Hold = Hold_sum / (double)(Hold_cnt);
		}
		if (Hnew_cnt != 0) {
			Hnew = Hnew_sum / (double)(Hnew_cnt);
		}
		bool res = ((Hold >= Hnew) || (dist01(rng) <= exp(Hold-Hnew)));
#if VERBOSE == 1
		//std::cout << i_old << " " << i_new << " " << Hold << " " << Hnew << " " << res << std::endl;
#endif
		return res;
	}

	void check(const size_t &step, const size_t &run) {}

	bool has_own_statistics() {
		return false;
	}
};

/**
 * TransitionMatrixSampler
 *
 * Following wang.j.02.transition.245 and shell.m.03.improved.9406.
 */
class TransitionMatrixSampler : public MCSampler {
private:
	const matrix_int_t& Q;
	vector_int_t column_sum;
public:
	TransitionMatrixSampler(boost::mt19937 &rng, const matrix_int_t& qmat, const settings_t &settings)
		: MCSampler(rng, settings), Q(qmat), column_sum(qmat.size1(),0) {
	}
	bool operator()(const int &E_old, const int &E_new) {
		column_sum[E_old]++;
#ifndef NDEBUG
		//int Hold_sum(0), Hnew_sum(0);
		//for (size_t i = 0; i < Q.size1(); i++) {
			//Hold_sum += Q(E_old,i);
			//Hnew_sum += Q(E_new,i);
		//}
		//assert(column_sum[E_old]==Hold_sum);
		//assert(column_sum[E_new]==Hnew_sum);
#endif
		if (column_sum[E_new] == 0 || column_sum[E_old] == 0) {
			return true;
		}
		double Tji, Tij;
		Tji = (double)Q(E_new, E_old) / column_sum[E_new];
		Tij = (double)Q(E_old, E_new) / column_sum[E_old];
		bool res = ((Tji >= Tij) || (dist01(rng) <= (Tji / Tij)));
		//std::cout << __LINE__ << " " << res << " " << Tji << " " << Tij << std::endl;
		return res;
	}

	void check(const size_t &step, const size_t &run) {}

	bool has_own_statistics() {
		return false;
	}
};


#endif /* end of include guard: MC_SAMPLER_H */
