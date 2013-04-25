#ifndef MC_SAMPLER_H

#define MC_SAMPLER_H


// Boost Random
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_01.hpp>

#include "typedefs.h"

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
};

class BoltzmannSampler : public MCSampler {
private:
	double kB;
	double T;
public:
	BoltzmannSampler(boost::mt19937 &rng, const matrix_int_t&, const settings_t &settings)
		: MCSampler(rng, settings) {
		kB = boost::any_cast<double>(settings.find("kB")->second);
		T  = boost::any_cast<double>(settings.find("temperature")->second);
	}

	bool operator()(const int &E_old, const int &E_new) {
		return ((E_old >= E_new) || (dist01(rng) <= exp((E_old-E_new)/(kB * T))));
	}

	bool operator()(const double &E_old, const double &E_new, const int &i_old, const int &i_new) {
		return ((E_old >= E_new) || (dist01(rng) <= exp((E_old-E_new)/(kB * T))));
	}

	void check(const size_t & step, const size_t &run) {
		//T = 2 - (1.9 * step/1e7);
		//T = sin(step/2e5)+1.1;
	}
};

#if 0
class WangLandauSampler : public MCSampler {
private:
	vector_int_t H;
	vector_double_t g;
	double ln_f;
	double flatness;
public:
	WangLandauSampler(boost::mt19937 &rng, const matrix_int_t&)
		: MCSampler(rng), H(ms), g(ms), ln_f(1.0), flatness(f) {
		H *= 0;
		g *= 0;
		for (size_t i = 0; i < ms; i++) {
			assert(H[i] == 0);
			assert(g[i] == 0);
		}
	}

	bool operator()(const int &E_old, const int &E_new) {
		using namespace std;
		bool res = ((g[E_old] >= g[E_new]) || (dist01(rng) <= exp(g[E_old]-g[E_new])));
		//cout<< __LINE__
		//	<< " (" << E_old << "," << E_new << ") "
		//	<< " [" << setw(12) << right << g[E_old] << "," << setw(12) << right << g[E_new] << "] "
		//	<< setw(12) << right << (g[E_old]-g[E_new])
		//	<< (res ? " T" : " F")
		//	<< setw(14) << right << ln_f
		//	<< g << " " << H
		//	<< std::endl;
		if (res) {
			H[E_new]++;
			g[E_new]+=ln_f;
		} else {
			H[E_old]++;
			g[E_old]+=ln_f;
		}
		return res;
	}

	void check(const size_t & step, const size_t &run) {
		int min = *std::min_element(H.begin(), H.end());
		if (min > 0 && (min > (0.99 * sum(H))/macro_states)) {
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
		//return true;
		return false;
	}

	//double calculate_error(const vector_double_t &exact) {
		//return ::calculate_error(exact, g, true);
	//}
};

class QualityMeasureBSampler : public MCSampler {
private:
	const matrix_int_t& Q;
public:
	QualityMeasureBSampler(boost::mt19937 &rng, size_t ms, double, double, double, const matrix_int_t& qmat)
		: MCSampler(rng, ms), Q(qmat) {
	}

	bool operator()(const int &E_old, const int &E_new) {
		size_t Hold_sum(0), Hnew_sum(0), Hold_cnt(0), Hnew_cnt(0);
		for (size_t i = 0; i < Q.size1(); i++) {
			if (Q(E_old,i) != 0) {
				Hold_sum += Q(E_old,i);
				Hold_cnt += 1;
			}
			if (Q(E_new,i) != 0) {
				Hnew_sum += Q(E_new,i);
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
		//std::cout << E_old << " " << E_new << " " << Hold << " " << Hnew << " " << res << std::endl;
#endif
		return res;
	}

	void check(const size_t &step, const size_t &run) {}

	bool has_own_statistics() {
		return false;
	}
};


class QualityMeasureASampler : public MCSampler {
private:
	const matrix_int_t& Q;
public:
	QualityMeasureASampler(boost::mt19937 &rng, size_t ms, double, double, double, const matrix_int_t& qmat)
		: MCSampler(rng, ms), Q(qmat) {
	}

	bool operator()(const int &E_old, const int &E_new) {
		size_t Hold_sum(0), Hnew_sum(0);
		for (size_t i = 0; i < Q.size1(); i++) {
			Hold_sum += Q(E_old,i);
			Hnew_sum += Q(E_new,i);
		}
		bool res = ((Hold_sum >= Hnew_sum) || (dist01(rng) <= exp(Hold_sum-Hnew_sum)));
#if VERBOSE == 1
		//std::cout << E_old << " " << E_new << " " << Hold << " " << Hnew << " " << res << std::endl;
#endif
		return res;
	}

	void check(const size_t &step, const size_t &run) {}

	bool has_own_statistics() {
		return false;
	}
};

#endif

#endif /* end of include guard: MC_SAMPLER_H */
