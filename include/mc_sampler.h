#ifndef MC_SAMPLER_H

#define MC_SAMPLER_H


// Boost Random
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_01.hpp>

#include "typedefs.h"
#include "q_matrix_tools.h"

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

	void check(const size_t & /*step*/, const size_t & /*run*/) {}

	bool has_own_statistics() {
		return false;
	}

	double calculate_error(const vector_double_t & /*exact*/) {
		return 0.0;
	}

	template<class T2>
	void rejected(const T2 &) {}
};


/**
 * BoltzmannSampler
 */
class BoltzmannSampler : public MCSampler {
private:
	double kB;
	double T;
	double beta;
public:
	BoltzmannSampler(boost::mt19937 &rng, const matrix_int_t&, const settings_t &settings)
		: MCSampler(rng, settings) {
		kB = boost::any_cast<double>(settings.find("kB")->second);
		T  = boost::any_cast<double>(settings.find("temperature")->second);
		beta = -1./(kB * T);
	}

	template<class T1>
	bool operator()(const T1 &dE) {
		return ((dE <= 0) || (dist01(rng) <= exp(beta * dE)));
	}

	template<class T1, class T2>
	bool operator()(const T1 &dE, const T2 & /*i_old*/, const T2 & /*i_new*/) {
		return ((dE <= 0) || (dist01(rng) <= exp(beta * dE)));
	}

	void check(const size_t & step, const size_t & /*run*/) {
		//T = 2 - (1.9 * step/1e7);
		//T = sin(step/2e5)+1.1;
	}
};

/**
 * WangLandauSampler
 */
class WangLandauSampler : public MCSampler {
protected:
	vector_int_t H;
	vector_double_t g;
	double ln_f;
	double flatness;
	size_t macrostates;
	size_t last_refinement_step, last_checked;
public:
	WangLandauSampler(boost::mt19937 &rng, const matrix_int_t&, const settings_t &settings)
		: MCSampler(rng, settings), ln_f(1.0), last_refinement_step(0), last_checked(0)
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
	bool operator()(const T1 & /*dE*/, const T2 &i_old, const T2 &i_new) {
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

	void check(const size_t & step, const size_t & /*run*/) {
		if ((step-last_checked) < 10*H.size() || (step-last_refinement_step) < 10*H.size()) {
			return;
		}
		last_checked = step;
		int minH = INT_MAX;
		int countH = 0;
		int sumH = 0;
		for (size_t i = 0; i < H.size(); ++i) {
			if (H[i] > 0) {
				countH++;
				sumH += H[i];
				if (H[i] < minH) {
					minH = H[i];
				}
			}
		}
		//std::cout << minH << " " <<  ((flatness * sumH)/countH) << " " <<  (minH > (flatness * sumH)/countH) << std::endl;
		if (minH > (flatness * sumH)/countH) {
			last_refinement_step = step;
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
		return rhab::calculate_error(exact, g, true);
	}
};

/**
 * QualityMeasureASampler
 */
class QualityMeasureASampler : public MCSampler {
private:
	boost::numeric::ublas::vector<int64_t> colsum;
public:
	QualityMeasureASampler(boost::mt19937 &rng, const matrix_int_t& qmat, const settings_t &settings)
		: MCSampler(rng, settings), colsum(qmat.size1()) {
		colsum *= 0;
	}

	template<class T1, class T2>
	bool operator()(const T1 & /*dE*/, const T2 &i_old, const T2 &i_new) {
		colsum[i_old]++;

		bool res = ((colsum[i_old] >= colsum[i_new]) || (dist01(rng) <= exp(colsum[i_old] - colsum[i_new])));

#if VERBOSE == 1
		int64_t Hold_sum(0), Hnew_sum(0);
		for (size_t i = 0; i < Q.size1(); i++) {
			Hold_sum += Q(i_old,i);
			Hnew_sum += Q(i_new,i);
		}
		assert(Hold_sum == colsum[i_old]);
		assert(Hnew_sum == colsum[i_new]);
		std::cout << i_old << " " << i_new << " " << Hold_sum << " " << Hnew_sum << " " << colsum[i_old] << " " << colsum[i_new] << " " << res << std::endl;
#endif
		return res;
	}

	template<class T2>
	void rejected(const T2 &i_old) {
		colsum[i_old]++;
	}

	void check(const size_t & /*step*/, const size_t & /*run*/) {}

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
	bool operator()(const T1 & /*dE*/, const T2 &i_old, const T2 &i_new) {
		int64_t Hold_sum(0), Hnew_sum(0), Hold_cnt(0), Hnew_cnt(0);
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

	void check(const size_t& /*step*/, const size_t& /*run*/) {}

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
	boost::numeric::ublas::vector<int64_t> colsum;
public:
	TransitionMatrixSampler(boost::mt19937 &rng, const matrix_int_t& qmat, const settings_t &settings)
		: MCSampler(rng, settings), Q(qmat), colsum(qmat.size1()) {
		colsum *= 0;
	}

	template<class T1, class T2>
	bool operator()(const T1& /*dE*/, const T2 &i_old, const T2 &i_new) {
		colsum[i_old]++;

		if (colsum[i_new] == 0 || colsum[i_old] == 0) {
			return true;
		}

		double Tji = (double)Q(i_new, i_old) / colsum[i_new];
		double Tij = (double)Q(i_old, i_new) / colsum[i_old];

		bool res = ((Tji >= Tij) || (dist01(rng) <= (Tji / Tij)));
#if VERBOSE == 1
		int64_t Hold_sum(0), Hnew_sum(0);
		for (size_t i = 0; i < Q.size1(); i++) {
			Hold_sum += Q(i_old,i);
			Hnew_sum += Q(i_new,i);
		}
		assert(Hold_sum == colsum[i_old]);
		assert(Hnew_sum == colsum[i_new]);
		std::cout << i_old << " " << i_new << " " << Hold_sum << " " << Hnew_sum << " " << colsum[i_old] << " " << colsum[i_new] << " " << res << std::endl;
#endif
		return res;
	}

	template<class T2>
	void rejected(const T2 &i_old) {
		colsum[i_old]++;
	}

	void check(const size_t& /*step*/, const size_t& /*run*/) {}

	bool has_own_statistics() {
		return false;
	}
};


#endif /* end of include guard: MC_SAMPLER_H */
