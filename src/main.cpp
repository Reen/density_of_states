#include <algorithm>
#include <fstream>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <vector>
#include <boost/algorithm/minmax_element.hpp>
#include <boost/format.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/math/common_factor.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/program_options.hpp>
#include <boost/cstdint.hpp>
#include <boost/tuple/tuple.hpp>

#include <libgen.h>

#include "config.h"
#include "GitSHA1.h"

#ifdef HAVE_DARWIN
#include <mach-o/dyld.h>
#endif


#define VERBOSE 0

typedef boost::numeric::ublas::vector<int> vector_int_t;
typedef boost::numeric::ublas::vector<double> vector_double_t;
typedef boost::numeric::ublas::matrix<int64_t> matrix_int_t;
typedef boost::numeric::ublas::matrix<double> matrix_double_t;

int get_gcd(const vector_int_t &a) {
	int curgcd = a[0];
	for (size_t i = 1; i < a.size()-1; i++) {
		curgcd = boost::math::gcd(curgcd, a[i]);
		if (curgcd == 1) {
			break;
		}
	}
	return curgcd;
}

double calculate_error(const vector_double_t &exact, const vector_double_t &dos, bool normalize = false) {
	vector_double_t::const_iterator i1 = dos.begin();
	vector_double_t::const_iterator i2 = exact.begin();
	double sum = 0.0;
	if (normalize) {
		double norm = 0;
		for (; i1 != dos.end(); i1++) {
			norm += exp(*i1);
		}
		i1 = dos.begin();
		for (; i1 != dos.end(); i1++, i2++) {
			sum += fabs((exp(*i1)/norm - *i2) / (*i2));
		}
		//std::cout << std::setprecision(22) << norm << " " << sum << " " << dos << " " << norm << std::endl;
	} else {
		for (; i1 != dos.end(); i1++, i2++) {
			sum += fabs((*i1 - *i2) / (*i2));
		}
	}
	return sum;
}

void normalize(vector_double_t &vec) {
	vec /= sum(vec);
}

vector_double_t calculate_dos_gth(matrix_double_t inner_mat) {
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

matrix_double_t normalize_q(const matrix_int_t & Q) {
	using namespace boost::numeric::ublas;
	matrix_double_t Qd(Q);
	for (size_t i = 0; i < Qd.size1(); i++) {
		double s = sum(row(Qd,i));
		if (s == 0) continue;
		row(Qd,i) /= s;
	}
	return Qd;
}

double calculate_error_q(const vector_double_t &exact, const matrix_double_t &Qd) {
	vector_double_t dos(calculate_dos_gth(Qd));
	return calculate_error(exact, dos);
}

double calculate_error_q_matrix(const matrix_double_t &Qex, const matrix_double_t &Q) {
	double value(0.0);
	for (size_t i = 0; i < Qex.data().size(); i++) {
		double tmp = Qex.data()[i] - Q.data()[i];
		value += tmp * tmp;
	}
	return value;
}

std::vector<size_t> build_histogram(const std::vector<double> &data, size_t bins, double &min, double &max, double &bin_width) {
	std::vector<size_t> hist(bins, 0);
	typedef std::vector<double>::const_iterator iterator;
	std::pair< iterator, iterator > res = boost::minmax_element(data.begin(), data.end());
	min = *(res.first);
	max = *(res.second)+std::numeric_limits<double>::epsilon();
	double diff = (max-min);
	bin_width = diff/bins;
	for (iterator i = data.begin(); i != data.end(); ++i) {
		size_t index = (*i-min)/bin_width;
		if (index < 0 || index >= hist.size()) {
			std::cerr << *i << " " << index << " " << min << " " << max << std::endl;
			throw std::out_of_range("Error: index of histogram out of bounds.");
		}
		hist[index]++;
	}
	return hist;
}

class MCSampler {
protected:
	boost::mt19937 &rng;
	boost::uniform_01<> dist01;
	size_t macro_states;
public:
	MCSampler(boost::mt19937 &_rng, size_t ms)
		: rng(_rng), macro_states(ms) {}

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
	BoltzmannSampler(boost::mt19937 &rng, size_t ms, double _kB, double _T, double, const matrix_int_t&)
		: MCSampler(rng, ms), kB(_kB), T(_T) {
	}

	bool operator()(const int &E_old, const int &E_new) {
		return ((E_old >= E_new) || (dist01(rng) <= exp((E_old-E_new)/(kB * T))));
	}

	void check(const size_t & step, const size_t &run) {
		//T = 2 - (1.9 * step/1e7);
		//T = sin(step/2e5)+1.1;
	}
};

class WangLandauSampler : public MCSampler {
private:
	vector_int_t H;
	vector_double_t g;
	double ln_f;
	double flatness;
public:
	WangLandauSampler(boost::mt19937 &rng, size_t ms, double, double, double f, const matrix_int_t&)
		: MCSampler(rng, ms), H(ms), g(ms), ln_f(1.0), flatness(f) {
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
		return true;
	}

	double calculate_error(const vector_double_t &exact) {
		return ::calculate_error(exact, g, true);
	}
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



class Simulation {
	//types
public:
	typedef std::vector< boost::tuple<size_t, double, double, double, double, double> > error_acc_t;

private:
	// variables set by parse_arguments
	size_t macro_states;
	size_t steps;
	size_t runs;
	size_t connections;
	size_t sampler;
	double T;
	double flatness;
	std::string tag;

	// variable set in constructor
	size_t error_check_f;
	std::vector<std::string> sampler_string;

	// random number generator
	boost::mt19937 rng;
	size_t seed;
	bool seed_set;

	// simulation variables and structures

	// mean & variance accumulators
	size_t error_acc_size;
	error_acc_t error_acc;
	std::vector<double> final_errors;

	// exact density of states
	vector_int_t dos_exact;

	// number of configurations
	size_t num_config;

	// exact density of states normalized
	vector_double_t dos_exact_norm;

	// map configuration to energy
	vector_int_t config_to_energy;

	// connectivity of states
	matrix_int_t mt;

	// microstate transition matrix
	matrix_double_t microstate_tm;

	// output stream
	std::ofstream out;

	// path to graph utility
	std::string graph_bin;

	// transition counting matrix
	matrix_int_t Q;

	// normalized Q matrix
	matrix_double_t Qd;

	// lumping matrix
	matrix_int_t lumper;

	// exact Q matrix
	matrix_double_t exact_q;

	bool parse_arguments(int argc, const char *argv[])
	{
		namespace po = boost::program_options;
		// parse argc / argv
		po::options_description desc("Allowed options");
		desc.add_options()
			("help", "produce help message")
			("macrostates,M", po::value<size_t>(&macro_states)->default_value(4), "Number of macrostates")
			("steps,S",       po::value<size_t>(&steps)->default_value(1000000),  "Number of steps per simulation")
			("runs,R",        po::value<size_t>(&runs)->default_value(1000),      "Number of simulations")
			("connections,C", po::value<size_t>(&connections)->default_value(3),  "Number of per-microstate connections")
			("temperature,T", po::value<double>(&T)->default_value(2.0),          "Temperature")
			("tag",           po::value<std::string>(&tag),                       "Additional tag to append to files")
			("seed",          po::value<size_t>(),                                "Random seed")
			("sampler",       po::value<size_t>(&sampler)->default_value(1),      "Sampler to use: 0 – Boltzmann, 1 – Wang-Landau, 2 – Q-Matrix")
			("flatness,f",    po::value<double>(&flatness)->default_value(0.99),  "Flatness parameter of the Wang-Landau algorithm")
			;
		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		if (vm.count("help")) {
			std::cout << "density of states toy "
				<< VERSION_MAJOR << "." << VERSION_MINOR
				<< "-" << g_GIT_SHA1 << std::endl;
			std::cout << desc << std::endl;
			std::exit(EXIT_FAILURE);
		}

		if (vm.count("seed")) {
			seed = vm["seed"].as<size_t>();
			seed_set = true;
		}
		//else {
		//}
		//std::cout << "random seed: " << seed << std::endl;
	}

	void initialize_rng()
	{
		if (!seed_set) {
			FILE* devran = fopen("/dev/urandom", "rb");
			size_t read = fread(&seed, sizeof(size_t), 1, devran);
			if (read != 1) {
				throw std::runtime_error("Could not read from /dev/urandom");
			}
			fclose(devran);
		}
		rng.seed(seed);
	}

	void safety_check()
	{
		if (connections < 2) {
			throw std::runtime_error("Error: number of connections is lower than two");
		}

		if (connections > num_config) {
			throw std::runtime_error("Error: number of connections is greater than number of configurations");
		}

		if (!has_regular_graph(num_config, connections)) {
			throw std::runtime_error("Error: Number of connections per micrstate and number of configurations don't fit.");
		}
	}

	void setup_variables()
	{
		// open output file
		std::string format;
		switch(sampler) {
		case 1:
			format = "dos_%1%_%2%S_%3%R_%4%M_%5%C_%7$0.2ff%8%%9%.out";
			break;
		case 2:
		case 3:
			format = "dos_%1%_%2%S_%3%R_%4%M_%5%C%8%%9%.out";
			break;
		case 0:
		default:
			format = "dos_%1%_%2%S_%3%R_%4%M_%5%C_%6$0.2fT%8%%9%.out";
		}
		std::string buf = str( boost::format(format)
				% sampler_string[sampler]
				% steps
				% runs
				% macro_states
				% connections
				% T
				% flatness
				% (tag.size() > 0 ? "_" : "")
				% tag
			);
		out.open( buf.c_str() );
		if (!out.good()) {
			throw std::runtime_error("Error: could not open output file");
		}

		if (steps % error_check_f != 0) {
			throw std::runtime_error("Error: check-frequency and number of steps don't match" );
		}
		error_acc_size = 9 * (size_t)log10(steps / error_check_f) + 1;
		error_acc.resize(error_acc_size);

		// calculate an exact density of states
		dos_exact = get_exact_dos(macro_states);
		num_config = sum(dos_exact);
		dos_exact_norm = dos_exact / static_cast<double>(num_config);
		config_to_energy = get_energy_map(dos_exact, num_config);

		out << "# dos_exact: " << dos_exact << '\n';
		out << "# dos_exact_norm: " << dos_exact_norm << '\n';
		out << "# config_to_energy: " << config_to_energy << '\n';

		// construct microstate graph
		mt.resize(num_config, connections, false);
		if (connections + 1 == num_config) {
			// short circuit around gengraph in case number of connections is num_config-1
			// gengraph is not able to generate such graphs before the heat death of the universe
			// all nodes simply have connection to all other nodes
			for (size_t i = 0; i < mt.size1(); i++) {
				size_t k = 0;
				for (size_t j = 0; j < mt.size2(); j++) {
					if (k == i) {
						k++;
					}
					mt(i,j) = k++;
				}
			}
		} else {
			// create cmd with arguments to run graph utility
			char buf[1024];
			boost::uniform_int<> graph_seed_dist(0, INT_MAX);
			int graph_seed = graph_seed_dist(rng);
			snprintf(buf, 1024, "echo %lu %lu | %s -s %i",
					connections, num_config, graph_bin.c_str(), graph_seed);

			out << "# " << buf << std::endl;
			FILE* fd = popen(buf, "r");
			if (!fd) {
				perror("Problem with pipe");
				throw std::runtime_error("Error: Can not execute graph utility!");
			}
			int row = 0;
			int col = 0;
			while (fgets(buf, 1024, fd)) {
				out << "# " << buf;
				std::istringstream in(buf);
				in >> row;
				for (size_t i = 0; i < connections; i++) {
					in >> col;
					mt(row, i) = col;
				}
			}
			pclose(fd);
		}

		out << "# mt: " << mt << '\n';

		microstate_tm =  boost::numeric::ublas::zero_matrix<matrix_int_t::value_type>(num_config, num_config);
		lumper = boost::numeric::ublas::zero_matrix<matrix_int_t::value_type>(num_config, macro_states);
		for (size_t i = 0; i < mt.size1(); i++) {
			out << "# ";
			for (size_t j = 0; j < mt.size2(); j++) {
				microstate_tm(i, mt(i,j)) = 1.0;//connections;
			}
			for (size_t j = 0; j < microstate_tm.size2(); j++) {
				out << std::setw(5) << microstate_tm(i,j);
			}
			out << '\n';
			//for (size_t j = 0; j < macro_states; j++) {
				//lumper(i, j) = (j == config_to_energy[i] ? 1 : 0);
			//}
			lumper(i, config_to_energy[i]) = 1;
		}

		matrix_double_t tmp = prod(trans(lumper),microstate_tm);
		exact_q = normalize_q(prod(tmp,lumper));
		out << "# microstate_tm: " << microstate_tm << '\n';
		out << "# lumper: " << lumper << '\n';
		out << "# exact_q: " << exact_q << std::endl;
	}

	vector_int_t get_exact_dos(size_t macro_states) {
		// use boost vectors ?
		int div;
		std::vector<double> dos_exact_d(macro_states);
		vector_int_t dos_exact_i(macro_states);
		for (size_t i = 1; i <= macro_states; i++) {
			dos_exact_d[i-1] = (i-(macro_states+1)/2.0);
			dos_exact_d[i-1] *= dos_exact_d[i-1];
			dos_exact_d[i-1] = -dos_exact_d[i-1] + (macro_states*macro_states+2*macro_states+1)/4.0;
			dos_exact_i[i-1] = dos_exact_d[i-1];
		}
		// compare _d and _i
		// divide by GCD until GCD is 1
		while((div = get_gcd(dos_exact_i)) != 1) {
			dos_exact_i /= div;
		}
		return dos_exact_i;
	}

	vector_int_t get_energy_map(const vector_int_t &dos_exact, const size_t &num_config) {
		vector_int_t energy_map(num_config);
		int k = 0;
		for (size_t i = 0; i < dos_exact.size(); i++) {
			for (int j = 0; j < dos_exact[i]; j++) {
				energy_map[k++] = i;
			}
		}
		return energy_map;
	}

	/**
	 * Determins if a regular graph can be constructed with the given parameters.
	 *
	 * @param size_t n number of nodes on the graph
	 * @param size_t r number of edges per node
	 */
	bool has_regular_graph(size_t n, size_t r) {
		// 3 nodes is minimum for a regular graph
		if (n < 3 || n <= r) return false;
		// we need at least 2 edges per node
		if (r == 2) return true;
		if (( r % 2 ) == 0) {
			return true;
		} else {
			// r is odd, so n must be even
			return (( n % 2 ) == 0);
		}
		return false;
	}

	template<typename Sampler>
	void mc_loop() {
		boost::uniform_int<> select_pos(0,connections-1);
		boost::uniform_01<> dist01;
		for (size_t run = 0; run < runs; run++) {
			error_check_f = 100;
			// start at random position
			size_t state = select_pos(rng);
			Q = boost::numeric::ublas::zero_matrix<int64_t>(macro_states, macro_states);
			//std::cout << "Q: " << Q << std::endl;
			Sampler sampler(rng, macro_states, kB, T, flatness, Q);
			size_t index = 0;
			for (size_t step = 1; step <= steps; step++) {
				size_t new_state = mt(state, select_pos(rng));
				int Eold = config_to_energy[state];
				int Enew = config_to_energy[new_state];
				Q(Eold, Enew)++;
				assert(new_state != state);
				if (sampler(Eold, Enew)) {
					state = new_state;
				}
				if (step > 0 && step % error_check_f == 0) {
					error_acc[index].get<0>() = step;
					Qd = normalize_q(Q);
					double err = calculate_error_q(dos_exact_norm, Qd);
					error_acc[index].get<1>() += err;
					error_acc[index].get<2>() += err*err;
					if (sampler.has_own_statistics()) {
						err = sampler.calculate_error(dos_exact_norm);
						error_acc[index].get<3>() += err;
						error_acc[index].get<4>() += err*err;
					}
					error_acc[index].get<5>() += calculate_error_q_matrix(exact_q, Qd);
					index++;
					if (step % (10*error_check_f) == 0) {
						error_check_f *= 10;
					}
				}
				sampler.check(step, run);
			}
			//std::cout << Q << std::endl;
			Qd = normalize_q(Q);
			//std::cout << Qd << std::endl;
			vector_double_t dos = calculate_dos_gth(Qd);
			final_errors.push_back(calculate_error(dos_exact_norm, dos));
			//std::cout << run << " " << dos << " " << calculate_error(dos_exact_norm, dos) << " " << sampler.calculate_error(dos_exact_norm) << std::endl;
		}
	}
public:
	static const double kB = 1.0;

	/**
	 *
	 */
	Simulation ()
		: error_check_f(100), seed_set(false)
	{
		// initialize list of available samplers
		sampler_string.push_back("BM");
		sampler_string.push_back("WL");
		sampler_string.push_back("QB");
		sampler_string.push_back("QA");

		// determin path to the "graph" executable
		// to add more OS see:
		// http://stackoverflow.com/questions/1023306/finding-current-executables-path-without-proc-self-exe
		// http://stackoverflow.com/questions/799679/programatically-retrieving-the-absolute-path-of-an-os-x-command-line-app
		char dir[1024] = ".";
#ifdef HAVE_LINUX
		ssize_t len = readlink("/proc/self/exe", dir, sizeof(dir)-1);
		if (len != -1) {
			dir[len] = '\0';
			dirname(dir); // dirname modifies argument
		} else {
			throw std::runtime_error("An error occured while determining the executable path");
		}
#endif
#ifdef HAVE_DARWIN
		char buf[1024];
		uint32_t size = sizeof(buf);
		if (_NSGetExecutablePath(buf, &size) == 0) {
			realpath(buf, dir);
			// dirname does not modify argument on OSX but instead returns a pointer
			// to internal memory. See man 3 dirname
			strcpy(dir, dirname(dir));
		} else {
			throw std::runtime_error("An error occured while determining the executable path");
		}
#endif
		graph_bin = std::string(dir) + "/graph";
#if VERBOSE==1
		std::cout << "Path to graph untility: " << graph_bin << std::endl;
#endif
	}

	int exec(int argc, const char *argv[]) {
		try {
			parse_arguments(argc, argv);
			initialize_rng();
			setup_variables();
			safety_check();
		} catch(std::exception& e) {
			std::cerr << "error: " << e.what() << std::endl;
			return EXIT_FAILURE;
		} catch(...) {
			std::cerr << "Exception of unknown type!" << std::endl;
			return EXIT_FAILURE;
		}

		// safety check
		if (steps % error_check_f != 0) {
			std::cerr << "Error: check-frequency and number of steps don't match" << std::endl;
			return EXIT_FAILURE;
		}

		if (run()) {
			return EXIT_SUCCESS;
		} else {
			return EXIT_FAILURE;
		}
	}

	bool run()
	{
		switch (sampler) {
		case 0:
			mc_loop<BoltzmannSampler>();
			break;
		case 1:
			mc_loop<WangLandauSampler>();
			break;
		case 2:
			mc_loop<QualityMeasureBSampler>();
			break;
		case 3:
			mc_loop<QualityMeasureASampler>();
			break;
		default:
			std::cerr << "Error: unknown sampler" << std::endl;
			return false;
		}

		out << "#\n#          time     mean_error      var_error   mean_error_s    var_error_s       qm_error\n";
		for (size_t i = 0; i < error_acc.size(); i++) {
			size_t time     = error_acc[i].get<0>();
			double mean     = error_acc[i].get<1>() / runs;
			double var_acc  = error_acc[i].get<2>() / runs;
			double mean2    = error_acc[i].get<3>() / runs;
			double var2_acc = error_acc[i].get<4>() / runs;
			double qm_err   = error_acc[i].get<5>() / runs;
			out << std::setw(15) << std::right
				<< time
				<< std::setw(15) << std::right
				<< mean
				<< std::setw(15) << std::right
				<< ((var_acc - mean * mean) / sqrt(runs))
				<< std::setw(15) << std::right
				<< mean2
				<< std::setw(15) << std::right
				<< ((var2_acc - mean2 * mean2) / sqrt(runs))
				<< std::setw(15) << std::right
				<< qm_err
				<< "\n";
		}

		out << "# Q:  " << Q  << "\n";
		out << "# Qd: " << Qd << "\n";
		out << "# diff (exact - Qd): "
			<< calculate_error_q_matrix(exact_q, Qd)
			<< "\n";

		double min, max, bin_width;
		std::vector<size_t> hist = build_histogram(final_errors, 20, min, max, bin_width);
		assert(final_errors.size() == std::accumulate(hist.begin(), hist.end(), 0));
		out << "#\n#\n# Histogram of errors:"
			<< "\n# min = " << min
			<< "\n# max = " << max
			<< "\n# bin_width = " << bin_width;
		for (size_t i = 0; i < hist.size(); i++) {
			out << "\n# "
				<< std::setw(12) << std::left << std::setprecision(6)
				<< (min + bin_width*i) << " | "
				<< std::setw(10) << hist[i] << " | ";
			for (size_t j = 0; j < (size_t)((double)hist[i]/final_errors.size()*100); j++) {
				out << '#';
			}
		}
		out << '\n';

	out.close();
		return true;
	}
};

int main(int argc, const char *argv[])
{
	Simulation s;
	return s.exec(argc, argv);
}

