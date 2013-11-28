#include "toydos_system.h"
#include "mc_sampler.h"
#include "rhab/misc.h"

#include <iomanip>

// Boost Bind
#include <boost/bind.hpp>

// Boost Format
#include <boost/format.hpp>

namespace po = boost::program_options;


/**
 * Private Methods:
 */

void ToyDosSystem::set_size(size_t L) {
	assert(L==macrostates);
	//macrostates = L;
	//lattice.resize(L, L);
	e_min = 1;
	e_max = L;
	n_bins = (e_max - e_min) + 1;
	assert(L==n_bins);
	Q.resize(n_bins, n_bins);
	Qd.resize(n_bins, n_bins);
}

void ToyDosSystem::set_gengraph_seed(int) {
	gengraph_seed_set = true;
}

void ToyDosSystem::setup_output() {
	std::vector<std::string> sampler_string;
	// initialize list of available samplers
	sampler_string.push_back("BM");
	sampler_string.push_back("WL");
	sampler_string.push_back("QB");
	sampler_string.push_back("QA");
	sampler_string.push_back("TM");

	std::string format;
	switch (sampler) {
		case 1:
			format = "toydos_%1%_%2%S_%3%R_%4%M_%5%C_%7$0.2ff%8%%9%.out";
			break;
		case 2:
		case 3:
		case 4:
			format = "toydos_%1%_%2%S_%3%R_%4%M_%5%C%8%%9%.out";
			break;
		case 0:
		default:
			format = "toydos_%1%_%2%S_%3%R_%4%M_%5%C_%6$0.2fT%8%%9%.out";
	}
	std::string buf = str( boost::format(format)
			% sampler_string[sampler]
			% steps
			% runs
			% n_bins
			% connections
			% boost::any_cast<double>(settings["temperature"])
			% boost::any_cast<double>(settings["flatness"])
			% (tag.size() > 0 ? "_" : "")
			% tag
		);
	out.open( buf.c_str() );
	if (!out.good()) {
		throw std::runtime_error("Error: could not open output file");
	}

}

void ToyDosSystem::setup_variables() {
	//out << "# Toy Dos Version: "
		//<< VERSION_MAJOR << "."
		//<< VERSION_MINOR << "-"
		//<< g_GIT_SHA1 << '\n';
	//out << "# seed: " << seed << '\n';

	// calculate an exact density of states
	dos_exact = get_exact_dos(macrostates);
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
		if (!gengraph_seed_set) {
			boost::uniform_int<> graph_seed_dist(0, INT_MAX);
			gengraph_seed = graph_seed_dist(rng.rng);
		}
		snprintf(buf, 1024, "echo %lu %lu | %s -s %i",
				connections, num_config, graph_bin.c_str(), gengraph_seed);

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

	out << "# mt: " << mt << std::endl;

	microstate_tm =  boost::numeric::ublas::zero_matrix<matrix_int_t::value_type>(num_config, num_config);
	lumper = boost::numeric::ublas::zero_matrix<matrix_int_t::value_type>(num_config, macrostates);
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
	exact_q.resize(macrostates,macrostates,false);
	rhab::normalize_q(prod(tmp,lumper), exact_q);
	out << "# microstate_tm: " << microstate_tm << '\n';
	out << "# lumper: " << lumper << '\n';
	out << "# exact_q: " << exact_q << std::endl;
}


void ToyDosSystem::safety_check() {

}


vector_int_t ToyDosSystem::get_exact_dos(size_t macro_states) {
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
	while((div = rhab::get_gcd(dos_exact_i)) != 1) {
		dos_exact_i /= div;
	}
	return dos_exact_i;
}

vector_int_t ToyDosSystem::get_energy_map(const vector_int_t &dos_exact, const size_t &num_config) {
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
bool ToyDosSystem::has_regular_graph(size_t n, size_t r) {
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



/**
 * Public Methods
 */

/**
 * Constructor
 */
ToyDosSystem::ToyDosSystem(settings_t &s) : SimulationSystem(s) {}

po::options_description ToyDosSystem::get_program_options() {
	po::options_description desc("ToyDos System Options");
	desc.add_options()
		(
			"macrostates,M",
			po::value<size_t>(&macrostates)->default_value(4)->notifier(boost::bind(&ToyDosSystem::set_size, this, _1)),
			"Number of macrostates"
		) (
			"connections,C",
			po::value<size_t>(&connections)->default_value(3),
			"Number of per-microstate connections"
		) (
			"gengraph_seed",
			po::value<int>(&gengraph_seed)->notifier(boost::bind(&ToyDosSystem::set_gengraph_seed, this, _1)),
			"Seed for gengraph tool"
		);
	return desc;
}

void ToyDosSystem::parse_arguments(boost::program_options::variables_map &vm) {
	if (vm.count("macrostates")) {
		settings["macrostates"] = vm["macrostates"].as<size_t>();
	}
}

void ToyDosSystem::setup() {
	SimulationSystem::setup();

	try {
		graph_bin     = boost::any_cast<std::string>(settings["executable_path"]) + "/graph";
		sampler       = boost::any_cast<size_t>(settings["sampler"]);
	} catch(const boost::bad_any_cast &e) {
		std::cerr << "Error while extracting settings in "
			<< __FILE__ << ":" << __LINE__ << "\nExecption: "
			<< e.what() << std::endl;
		throw e;
	}

	size_t error_acc_size = 9 * (size_t)log10(steps / error_check_f) + 1;
	error_acc.resize(error_acc_size);

	//read_exact_dos();
	setup_output();
	setup_variables();
	safety_check();
}

bool ToyDosSystem::run() {
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
	case 4:
		mc_loop<TransitionMatrixSampler>();
		break;
	default:
		std::cerr << "Error: unknown sampler" << std::endl;
		return false;
	}

	write_output();

	return true;
}
