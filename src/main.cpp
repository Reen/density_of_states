#include <algorithm>
#include <fstream>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <vector>
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
#include "config.h"
#include "GitSHA1.h"

#define VERBOSE 0

namespace po = boost::program_options;
namespace ublas = boost::numeric::ublas;
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

vector_int_t get_energy_map(const vector_int_t &dos_exact, const size_t &Nconfig) {
	vector_int_t energy_map(Nconfig);
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

void normalize(vector_double_t &vec) {
	vec /= sum(vec);
}

vector_double_t calculate_dos_gth(matrix_double_t inner_mat) {
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

double calculate_error_q(const vector_double_t &exact, const matrix_int_t &Q) {
	matrix_double_t Qd(normalize_q(Q));
	vector_double_t dos(calculate_dos_gth(Qd));
	return calculate_error(exact, dos);
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
	BoltzmannSampler(boost::mt19937 &rng, size_t ms, double _kB, double _T)
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
public:
	WangLandauSampler(boost::mt19937 &rng, size_t ms, double, double)
		: MCSampler(rng, ms), H(ms), g(ms), ln_f(1.0) {
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
		/*cout<< __LINE__
			<< " (" << E_old << "," << E_new << ") "
			<< " [" << setw(12) << right << g[E_old] << "," << setw(12) << right << g[E_new] << "] "
			<< setw(12) << right << (g[E_old]-g[E_new])
			<< (res ? " T" : " F")
			<< setw(14) << right << ln_f
			<< g << " " << H
			<< std::endl;*/
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

template<typename Sampler>
void mc_loop(boost::mt19937 &rng,
		const size_t &macro_states,
		const size_t &steps,
		const size_t &runs,
		const size_t &connections,
		const double &T,
		const double &kB,
		const size_t &error_check_f,
		vector_double_t &error_acc,
		vector_double_t &error_acc2,
		vector_double_t &error_s_acc,
		vector_double_t &error_s_acc2,
		vector_int_t &config_to_energy,
		matrix_int_t &mt,
		vector_double_t &dos_exact_norm
		) {
	boost::uniform_int<> select_pos(0,connections-1);
	boost::uniform_01<> dist01;
	for (size_t run = 0; run < runs; run++) {
		// start at random position
		size_t state = select_pos(rng);
		matrix_int_t Q(boost::numeric::ublas::zero_matrix<int64_t>(macro_states, macro_states));
		//std::cout << "Q: " << Q << std::endl;
		Sampler sampler(rng, macro_states, kB, T);
		for (size_t step = 0; step < steps; step++) {
			size_t new_state = mt(state, select_pos(rng));
			int Eold = config_to_energy[state];
			int Enew = config_to_energy[new_state];
			Q(Eold, Enew)++;
			assert(new_state != state);
			if (sampler(Eold, Enew)) {
				state = new_state;
			}
			if (step % error_check_f == 0) {
				double err = calculate_error_q(dos_exact_norm, Q);
				error_acc[step / error_check_f]  += err;
				error_acc2[step / error_check_f] += err*err;
				if (sampler.has_own_statistics()) {
					err = sampler.calculate_error(dos_exact_norm);
					error_s_acc[step / error_check_f]  += err;
					error_s_acc2[step / error_check_f] += err*err;
				}
			}
			sampler.check(step, run);
		}
		//std::cout << Q << std::endl;
		matrix_double_t Qd(normalize_q(Q));
		//std::cout << Qd << std::endl;
		vector_double_t dos = calculate_dos_gth(Qd);
		//std::cout << run << " " << dos << " " << calculate_error(dos_exact_norm, dos) << " " << sampler.calculate_error(dos_exact_norm) << std::endl;
	}
}

int main(int argc, const char *argv[])
{
	// Variables set by program_options
	size_t macro_states = 4;
	size_t steps = 10000000;
	size_t runs = 1000;
	size_t connections;
	double T = 2.0;
	std::string tag;
	const double kB = 1.0;
	const size_t error_check_f = 1000;
	size_t sampler = 0;
	boost::mt19937 rng;
	std::vector<std::string> sampler_string;
	sampler_string.push_back("BM");
	sampler_string.push_back("WL");

	// parse argc / argv
	try {
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
			("sampler",       po::value<size_t>(&sampler)->default_value(1),      "Sampler to use: 0 - Boltzmann, 1 - WangLandau")
			;
		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		if (vm.count("help")) {
			std::cout << "density of states toy "
				<< VERSION_MAJOR << "." << VERSION_MINOR
				<< "-" << g_GIT_SHA1 << std::endl;
			std::cout << desc << std::endl;
			return 1;
		}

		size_t seed;
		if (vm.count("seed")) {
			seed = vm["help"].as<size_t>();
		} else {
			FILE* devran = fopen("/dev/urandom", "rb");
			fread(&seed, sizeof(size_t), 1, devran);
			fclose(devran);
		}
		std::cout << "random seed: " << seed << std::endl;
		rng.seed(seed);
	} catch(std::exception& e) {
		std::cerr << "error: " << e.what() << "\n";
		return 1;
	} catch(...) {
		std::cerr << "Exception of unknown type!\n";
	}

	// safety check
	if (steps % error_check_f != 0) {
		std::cerr << "Error: check-frequency and number of steps don't match" << std::endl;
		exit(1);
	}

	// mean & variance accumulators
	vector_double_t error_acc(steps / error_check_f);
	vector_double_t error_acc2(steps / error_check_f);
	vector_double_t error_s_acc(steps / error_check_f);
	vector_double_t error_s_acc2(steps / error_check_f);

	// set to zero
	error_acc *= 0;
	error_acc2 *= 0;
	error_s_acc *= 0;
	error_s_acc2 *= 0;

	// calculate an exact density of states
	vector_int_t dos_exact = get_exact_dos(macro_states);
	size_t Nconfig = sum(dos_exact);
	vector_double_t dos_exact_norm = dos_exact / static_cast<double>(Nconfig);
	vector_int_t config_to_energy = get_energy_map(dos_exact, Nconfig);

	if (connections < 2) {
		std::cerr << "Error: number of connections is lower than two" << std::endl;
		return 1;
	}

	if (connections > Nconfig) {
		std::cerr << "Error: number of connections is greater than number of configurations" << std::endl;
		return 1;
	}

	if (!has_regular_graph(Nconfig, connections)) {
		std::cerr << "Error: Number of connections per micrstate and number of configurations don't fit." << std::endl;
		return 1;
	}

	// output file
	std::ofstream out;
	{
		char buf[1024];
		snprintf(buf, 1024, "dos_%s_%luS_%luR_%luM_%luC_%fT%s%s.out",
				sampler_string[sampler].c_str(),
				steps,
				runs,
				macro_states,
				connections,
				T,
				(tag.size() > 0 ? "_" : ""),
				tag.c_str()
				);
		out.open(buf);
	}
	out << "# dos_exact: " << dos_exact << std::endl;
	out << "# dos_exact_norm: " << dos_exact_norm << std::endl;
	out << "# config_to_energy: " << config_to_energy << std::endl;

	// construct microstate graph
	matrix_int_t mt(Nconfig, connections);
	if (connections + 1 == Nconfig) {
		// short circuit around gengraph in case number of connections is Nconfig-1
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
		char buf[1024];
		int graph_seed = rng();
		snprintf(buf, 1024, "echo %lu %lu | ./graph -s %i", connections, Nconfig, graph_seed);
		out << "# " << buf << std::endl;
		FILE* fd = popen(buf, "r");
		//FILE* fd = popen("ps aux", "r");
		if (!fd) {
			perror("Problem with pipe");
			return 1;
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

	switch (sampler) {
	case 0:
		mc_loop<BoltzmannSampler>(rng, macro_states, steps, runs, connections, T, kB,
			error_check_f, error_acc, error_acc2, error_s_acc, error_s_acc2,
			config_to_energy, mt, dos_exact_norm);
		break;
	case 1:
		mc_loop<WangLandauSampler>(rng, macro_states, steps, runs, connections, T, kB,
			error_check_f, error_acc, error_acc2, error_s_acc, error_s_acc2,
			config_to_energy, mt, dos_exact_norm);
		break;
	default:
		std::cerr << "Error: unknown sampler" << std::endl;
	}
	error_acc /= runs;
	error_acc2 /= runs;
	for (size_t i = 0; i < error_acc.size(); i++) {
		out << std::setw(15) << std::right << (i * error_check_f)
			<< std::setw(15) << std::right << error_acc[i]
			<< std::setw(15) << std::right
			<< ((error_acc2[i] - error_acc[i] * error_acc[i]) / sqrt(runs))
			<< std::setw(15) << std::right << error_s_acc[i]
			<< std::setw(15) << std::right
			<< ((error_s_acc2[i] - error_s_acc[i] * error_s_acc[i]) / sqrt(runs))
			<< "\n";
	}

	out.close();

	return 0;
}
