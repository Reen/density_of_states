#ifndef MC_SAMPLER_H

#define MC_SAMPLER_H

#include <iomanip>


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

	void get_parameter(double& /*param*/) const {}

	double calculate_error(const vector_double_t& /*exact*/,
						   error_mat_tuple_t /*error_matrices*/,
						   const size_t& /*index*/) {
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
	size_t steps;
	size_t schedule;
public:
	BoltzmannSampler(boost::mt19937 &rng, const matrix_int_t&, const settings_t &settings)
		: MCSampler(rng, settings) {
		kB    = boost::any_cast<double>(settings.find("kB")->second);
		T     = boost::any_cast<double>(settings.find("temperature")->second);
		steps = boost::any_cast<size_t>(settings.find("steps")->second);
		schedule = boost::any_cast<size_t>(settings.find("schedule")->second);
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
		double Tcur;
		switch(schedule) {
		// Linear Schedule
		case 1:
			Tcur = (0.1-T)/steps*step+T;
			beta = -1./(kB * Tcur);
			break;
		// Sine Schedule
		case 2:
			Tcur = T*(sin(step/(steps/100))+1.0);
			beta = -1./(kB * Tcur);
			break;
		case 0:
		default:
			break;
		}
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
	size_t last_refinement_step, last_checked, max_count0;
public:
	WangLandauSampler(boost::mt19937 &rng, const matrix_int_t&, const settings_t &settings)
		: MCSampler(rng, settings), ln_f(1.0), last_refinement_step(0), last_checked(0)
	{
		macrostates = boost::any_cast<size_t>(settings.find("macrostates")->second);
		flatness    = boost::any_cast<double>(settings.find("flatness")->second);
		max_count0  = boost::any_cast<size_t>(settings.find("num_unvisited_energies")->second);

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

	void check(const size_t & step, const size_t & run) {
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

		// We have not explored all available states
		// in this refinement step => return.
		if (H.size() - countH != max_count0) {
			return;
		}

		//std::cout << minH << " " <<  ((flatness * sumH)/countH) << " " <<  (minH > (flatness * sumH)/countH) << " " << ln_f << std::endl;
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
	}

	void get_parameter(double& param) const {
		param = ln_f;
	}

	double calculate_error(const vector_double_t& exact,
						   error_mat_tuple_t error_matrices,
						   const size_t& index) {
		return rhab::calculate_error(exact, g, error_matrices.get<3>(), index, true);
	}
};

class WangLandau1tSampler : public WangLandauSampler {
protected:
	size_t count0_pre1, count0_pre2, step_pre2;
	double param_c;
	bool use_one_t;
	size_t step_, param_s;
public:
	WangLandau1tSampler(boost::mt19937& rng, const matrix_int_t& Q, const settings_t& settings)
		: WangLandauSampler(rng, Q, settings), count0_pre1(0), count0_pre2(0),
		  step_pre2(0), use_one_t(false) {
		param_c = boost::any_cast<double>(settings.find("one-over-t-c")->second);
		param_s = boost::any_cast<size_t>(settings.find("one-over-t-s")->second);
	}

	template<class T1, class T2>
	bool operator()(const T1 & /*dE*/, const T2 &i_old, const T2 &i_new) {
		bool res = ((g[i_old] >= g[i_new]) || (dist01(rng) <= exp(g[i_old]-g[i_new])));

		if (use_one_t) {
			// step / H.size() is the "time" \f$t\f$ as defined by Belardinelli et al.
			ln_f = param_c/( step_/static_cast<double>(H.size()) );
		}

		if (res) {
			H[i_new]++;
			g[i_new]+=ln_f;
		} else {
			H[i_old]++;
			g[i_old]+=ln_f;
		}
		return res;
	}

	void check(const size_t & step, const size_t & run) {
		step_ = step;
		if ((step-last_checked) < 100*H.size()
			|| (step-last_refinement_step) < 100*H.size()
			|| use_one_t) {
			return;
		}
		last_checked = step;
		size_t count0 = 0;
		//int minH = INT_MAX;
		//int countH = 0;
		//int sumH = 0;
		for (size_t i = 0; i < H.size(); ++i) {
			if (H[i] == 0) {
				count0++;
			}
		}

		if (count0 == max_count0
			&& count0 == count0_pre2) {
			last_refinement_step = step;
			H *= 0;
			ln_f /= 2.0;

			// step / H.size() is the "time" \f$t\f$ as defined by Belardinelli et al.
			if (ln_f <= param_c/( step/static_cast<double>(H.size()) )) {
				use_one_t = true;
				//std::cout << "switching to 1/t after step " << step << std::endl;
			}

#if VERBOSE == 1
			std::cout << std::setw(15) << std::right << run
					  << std::setw(15) << std::right << step
					  << std::setw(15) << std::right << ln_f
					  << std::setw(15) << std::right << count0
					  << std::setw(15) << std::right << count0_pre1
					  << std::setw(15) << std::right << count0_pre2
					  << "\n";
#endif
		}
		if (step-step_pre2 > param_s * H.size()) {
			count0_pre2 = count0_pre1;
			step_pre2 = step;
		}
		count0_pre1 = count0;
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
			std::fill(colsum.begin(), colsum.end(), 0);
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
		std::fill(colsum.begin(), colsum.end(), 0);
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
