#include "toydos_system.h"
#include "mc_sampler.h"
#include "rhab/misc.h"

#include <iomanip>
#include <set>

// Boost Bind
#include <boost/bind.hpp>

// Boost Format
#include <boost/format.hpp>

// Boost Random Uniform_smallint
#include <boost/random/uniform_smallint.hpp>

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

  std::fill(Q.data().begin(), Q.data().end(), 0);
  std::fill(Qd.data().begin(), Qd.data().end(), 0.0);
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
  sampler_string.push_back("1t");

  std::string format;
  switch (sampler) {
    case 1:
      format = "toydos_%1%_%2%S_%3%R_%4%M_%5%C_%11%-%12%L_%7$0.2ff%8%%9%.out";
      break;
    case 2:
    case 3:
    case 4:
    case 5:
      format = "toydos_%1%_%2%S_%3%R_%4%M_%5%C_%11%-%12%L%8%%9%.out";
      break;
    case 0:
      format = "toydos_%1%_%2%S_%3%R_%4%M_%5%C_%11%-%12%L_%6$0.2fT_%10%Sd%8%%9%.out";
      break;
    default:
      throw std::runtime_error("Unknown smapler");
  }
  boost::format fmt(format);
  fmt.exceptions(boost::io::all_error_bits ^ (
        boost::io::too_many_args_bit | boost::io::too_few_args_bit
        ));
  fn_template = str( fmt
            /*  1 */ % sampler_string[sampler]
            /*  2 */ % steps
            /*  3 */ % "%1%"
            /*  4 */ % n_bins
            /*  5 */ % connections
            /*  6 */ % boost::any_cast<double>(settings["temperature"])
            /*  7 */ % boost::any_cast<double>(settings["flatness"])
            /*  8 */ % (tag.size() > 0 ? "_" : "")
            /*  9 */ % tag
            /* 10 */ % boost::any_cast<size_t>(settings["schedule"])
            /* 11 */ % num_layers
            /* 12 */ % num_connections_between_layers
      );

  out << "#";
  out << "\n# sampler:        " << sampler << " - " << sampler_string[sampler];
  out << "\n# steps:          " << steps;
  out << "\n# n_bins:         " << n_bins;
  out << "\n# macrostates:    " << macrostates;
  out << "\n# connections:    " << connections;
  out << "\n# temperature:    " << boost::any_cast<double>(settings["temperature"]);
  out << "\n# flatness:       " << boost::any_cast<double>(settings["flatness"]);
  out << "\n# one-over-t-c:   " << boost::any_cast<double>(settings["one-over-t-c"]);
  out << "\n# one-over-t-s:   " << boost::any_cast<size_t>(settings["one-over-t-s"]);
  if (gengraph_seed_set) {
    out << "\n# grengraph_seed: " << gengraph_seed;
  }
  out << "\n# reduce_micro:   " << reduce_microstates;
  out << "\n# layers:         " << num_layers;
  out << "\n# layer_conn:     " << num_connections_between_layers;
  out << "\n";
}


matrix_int_t ToyDosSystem::generate_single_graph(const size_t& nconfig,
                                                 const size_t& conn,
                                                 const int& gengraph_seed_offset) {
  matrix_int_t mt_tmp(nconfig, conn);
  if (conn + 1 == nconfig) {
    // Short circuit around \verb!gengraph! in case the number of connections is
    // nconfig-1. \verb!gengraph! is then not able to generate such graphs
    // before the heat death of the universe.
    // All nodes simply have connection to all other nodes.
    for (size_t i = 0; i < mt_tmp.size1(); i++) {
      size_t k = 0;
      for (size_t j = 0; j < mt_tmp.size2(); j++) {
        if (k == i) {
          k++;
        }
        mt_tmp(i,j) = k++;
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
        conn, nconfig, graph_bin.c_str(),
        gengraph_seed + gengraph_seed_offset);

    out << "# " << buf << std::endl;
    FILE* fd = popen(buf, "r");
    if (!fd) {
      perror("Problem with pipe");
#ifdef USE_MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
#endif
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
#ifndef NDEBUG
        if (col < 0 || col > (int)(nconfig)-1) {
          std::cerr << buf << std::flush;
          throw(std::runtime_error("Error: Graph utility gave unexpected results"));
        }
#endif
        mt_tmp(row, i) = col + (gengraph_seed_offset * nconfig);
      }

      // I like to keep them in order
      std::sort(&mt_tmp(row, 0), &mt_tmp(row, connections-1));
    }
    if ((size_t)(row+1) != mt_tmp.size1()) {
#ifdef USE_MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
#endif
      throw std::runtime_error("Error reading graph utility results!");
    }
    pclose(fd);
  }
  return mt_tmp;
}

matrix_int_t ToyDosSystem::generate_multi_layer_graph() {
  if (num_config % num_layers != 0) {
    throw std::runtime_error("Error: number of microstates is not divisible \
        by the number of layers. Try to set --reduce_microstates to false.");
  }
  matrix_int_t mt_tmp(num_config, connections);
  size_t num_config_per_layer = num_config/num_layers;
  for (size_t i = 0; i < num_layers; ++i) {
    matrix_int_t mt_layer = generate_single_graph(num_config_per_layer,
                                                  connections, i);
    for (size_t j = 0; j < mt_layer.size1(); j++) {
      for (size_t k = 0; k < mt_layer.size2(); ++k) {
        mt_tmp(j + (num_config_per_layer*i), k) = mt_layer(j, k);
      }
    }
  }

  boost::uniform_smallint<> state_dist(0, num_config_per_layer-1);
  boost::uniform_smallint<> conn_dist(0, connections-1);

  // combine
  std::set<int> connected;
  for (size_t i = 0; i < num_layers-1; ++i) {
    // connect layer i to i+1

    for (size_t j = 0; j < num_connections_between_layers; ++j) {
      int c_fwd_i;
      // select a microstate to start from that has
      // not been connected to another layer yet
      do {
        c_fwd_i = state_dist(rng.rng)
          + (num_config_per_layer * i);
        // removing the previous line leads to a star-like structure
        // for the set of clusters
      } while (connected.count(c_fwd_i) > 0);
      connected.insert(c_fwd_i);

      const int c_fwd_j_idx   = conn_dist(rng.rng);
      const int c_fwd_j       = mt_tmp(c_fwd_i, c_fwd_j_idx);

      const int c_bwd_i       = c_fwd_j;
      const int c_bwd_j_idx   = std::find(&mt_tmp(c_bwd_i, 0),
                                          &mt_tmp(c_bwd_i, connections-1),
                                          c_fwd_i) - &mt_tmp(c_bwd_i, 0);

      int c_fwd_ip1;
      // select a microstate to start from that has
      // not been connected to another layer yet
      do {
        c_fwd_ip1 = state_dist(rng.rng) + (num_config_per_layer * (i+1));
      } while (connected.count(c_fwd_ip1) > 0);
      connected.insert(c_fwd_ip1);

      const int c_fwd_jp1_idx = conn_dist(rng.rng);
      const int c_fwd_jp1     = mt_tmp(c_fwd_ip1, c_fwd_jp1_idx);

      const int c_bwd_ip1     = c_fwd_jp1;
      const int c_bwd_jp1_idx = std::find(&mt_tmp(c_bwd_ip1, 0),
          &mt_tmp(c_bwd_ip1, connections-1),
          c_fwd_ip1) - &mt_tmp(c_bwd_ip1, 0);

      std::swap(mt_tmp(c_fwd_i, c_fwd_j_idx), mt_tmp(c_fwd_ip1, c_fwd_jp1_idx));
      std::swap(mt_tmp(c_bwd_i, c_bwd_j_idx), mt_tmp(c_bwd_ip1, c_bwd_jp1_idx));
    }
  }
  return mt_tmp;
}

void ToyDosSystem::setup_variables() {
  // calculate an exact density of states
  dos_exact = get_exact_dos(macrostates);
  num_config = sum(dos_exact);
  dos_exact_norm = dos_exact / static_cast<double>(num_config);
  config_to_energy = get_energy_map(dos_exact, num_config);

  out << "# dos_exact: " << dos_exact << '\n';
  out << "# dos_exact_norm: " << dos_exact_norm << '\n';
  out << "# config_to_energy: " << config_to_energy << '\n';

  // @todo: check that num_config is divisible by num_layers
  // construct microstate graph
  if (world_rank == 0) {
    if (num_layers == 1) {
      mt = generate_single_graph(num_config, connections, 0);
    } else {
      mt = generate_multi_layer_graph();
    }
  }
#ifdef USE_MPI
    size_t mt_size_1 = mt.size1();
    size_t mt_size_2 = mt.size2();
    MPI_Bcast(&mt_size_1, (sizeof(size_t)), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mt_size_2, (sizeof(size_t)), MPI_BYTE, 0, MPI_COMM_WORLD);
    if (world_rank != 0) {
      mt.resize(mt_size_1, mt_size_2, false);
    }
    MPI_Bcast(&mt.data()[0], (sizeof(matrix_int_t::value_type))*mt_size_1*mt_size_2,
              MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

  out << "# mt: " << mt << std::endl;
  if (print_microstate_graph) {
    matrix_int_t microstate_tm(num_config, num_config);
    microstate_tm.clear();
    for (size_t i = 0; i < mt.size1(); i++) {
      for (size_t j = 0; j < mt.size2(); j++) {
        microstate_tm(i, mt(i,j)) = 1;//connections;
      }
      for (size_t j = 0; j < microstate_tm.size2(); j++) {
        out << std::setw(2) << microstate_tm(i,j);
      }
      out << '\n';
    }
  }
  exact_q.resize(macrostates, macrostates, false);
  /*
  // Variant 1: (actually produces lumper and microstate_tm
  cmatrix_int_t microstate_tm(num_config, num_config);
  cmatrix_int_t lumper(num_config, macrostates);
  exact_q.resize(macrostates,macrostates,false);
  for (size_t i = 0; i < mt.size1(); i++) {
  for (size_t j = 0; j < mt.size2(); j++) {
  microstate_tm(i, mt(i,j)) = 1;//connections;
  }
  if (num_config < 200) {
  out << "# ";
  for (size_t j = 0; j < microstate_tm.size2(); j++) {
  out << std::setw(5) << microstate_tm(i,j);
  }
  out << '\n';
  }
  lumper(i, config_to_energy[i]) = 1;
  }

  cmatrix_int_t tmp = prod(trans(lumper), microstate_tm);
  rhab::normalize_q(prod(tmp, lumper), exact_q);
  if (num_config < 200) {
  out << "# microstate_tm: " << microstate_tm << '\n';
  out << "# lumper: " << lumper << '\n';
  }
  out << "# exact_q: " << exact_q << std::endl;
  */

  // variant 2:
  matrix_int_t qq(macrostates, macrostates);
  qq.clear();
  for (size_t i = 0; i < mt.size1(); i++) {
    for (size_t j = 0; j < mt.size2(); j++) {
      qq(config_to_energy[i], config_to_energy[mt(i, j)])++;
    }
  }
  exact_q.clear();
  rhab::normalize_q(qq, exact_q);

  out << "# exact_q: " << exact_q << std::endl;
}


void ToyDosSystem::safety_check() {

}


vector_int_t ToyDosSystem::get_exact_dos(size_t macro_states) {
  // use boost vectors ?
  //std::vector<double> dos_exact_d(macro_states);
  vector_int_t dos_exact_i(macro_states);
  /*for (size_t i = 1; i <= macro_states; i++) {
    dos_exact_d[i-1] = (i-(macro_states+1)/2.0);
    dos_exact_d[i-1] *= dos_exact_d[i-1];
    dos_exact_d[i-1] = -dos_exact_d[i-1] + (macro_states*macro_states+2*macro_states+1)/4.0;
    dos_exact_i[i-1] = dos_exact_d[i-1];
  }*/
  // compare _d and _i
  // divide by GCD until GCD is 1
  /*if (reduce_microstates) {
    int div;
    while((div = rhab::get_gcd(dos_exact_i)) != 1) {
      dos_exact_i /= div;
    }
  }*/
  for (size_t i = 0; i < macro_states; i++) {
    dos_exact_i[i] = floor(pow(3., 1.+pow(100.*i, 0.305)));
    //std::cout << i << " " << dos_exact_i[i] << std::endl;
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
ToyDosSystem::ToyDosSystem(settings_t &s) : SimulationSystem(s), gengraph_seed_set(false) {}

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
        ) (
          "reduce_microstates",
          po::value<bool>(&reduce_microstates)->default_value(true),
          "Whether to reduce the number of microstates by ggT. Should be set to 0 for layers > 1."
          ) (
            "print_microstate_graph",
            po::value<bool>(&print_microstate_graph)->default_value(false),
            "Whether to print the microstate graph to cout"
            ) (
              "layers",
              po::value<size_t>(&num_layers)->default_value(1),
              "Number of layers into which the microstates are clustered."
              ) (
                "layer_connections",
                po::value<size_t>(&num_connections_between_layers)->default_value(1),
                "Number of connections between layers"
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

  size_t num_unvisited_energies = 0;
  settings["num_unvisited_energies"] = num_unvisited_energies;

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
    case 5:
      mc_loop<WangLandau1tSampler>();
      break;
    default:
      std::cerr << "Error: unknown sampler" << std::endl;
      return false;
  }

  return true;
}

/* vim: set ts=2 sw=2 sts=2 tw=0 expandtab :*/
