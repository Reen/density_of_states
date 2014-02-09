#include "lennard_jones_system.h"
#include "mc_sampler.h"

// C++ Standard Library
#include <iomanip>

// Boost Bind
#include <boost/bind.hpp>

// Boost Format
#include <boost/format.hpp>

namespace po = boost::program_options;


/**
 * Private Methods:
 */

void LennardJonesSystem::set_bins(size_t nb) {
  n_bins = nb;
  Q.resize(n_bins, n_bins);
  Qd.resize(n_bins, n_bins);

  std::fill(Q.data().begin(), Q.data().end(), 0);
  std::fill(Qd.data().begin(), Qd.data().end(), 0.0);
}

void LennardJonesSystem::set_particles(size_t N) {
  particles.resize(N);
  for (size_t i = 0; i < particles.size(); i++) {
    for (size_t j = 0; j < 3; j++) {
      particles[i][j] = size * dist01(rng.rng);
    }
  }

  ran_particle.param(boost::uniform_smallint<size_t>::param_type(0,N-1));
}

double calculate_distance(const boost::array<double, 3> &a, const boost::array<double, 3> &b) {
  double r = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    double tmp = a[i] - b[i];
    r += tmp*tmp;
  }
  return sqrt(r);
}

void LennardJonesSystem::setup_output() {
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
      format = "lj_%1%_%2%S_%3%R_%4%M_%6$0.2ff%7%%8%.out";
      break;
    case 2:
    case 3:
    case 4:
    case 5:
      format = "lj_%1%_%2%S_%3%R_%4%M%7%%8%.out";
      break;
    case 0:
      format = "lj_%1%_%2%S_%3%R_%4%M_%5$0.2fT_%9%Sd%7%%8%.out";
      break;
    default:
      throw std::runtime_error("Unknown smapler");
  }
  boost::format fmt(format);
  fmt.exceptions(boost::io::all_error_bits ^ (
        boost::io::too_many_args_bit | boost::io::too_few_args_bit
        ));
  std::string buf = str( fmt
      % sampler_string[sampler]
      % steps
      % runs
      % n_bins
      % boost::any_cast<double>(settings["temperature"])
      % boost::any_cast<double>(settings["flatness"])
      % (tag.size() > 0 ? "_" : "")
      % tag
      % boost::any_cast<size_t>(settings["schedule"])
    );

  open_output_files(buf);

  write_header();

  out << "#";
  out << "\n# sampler:      " << sampler << " - " << sampler_string[sampler];
  out << "\n# steps:        " << steps;
  out << "\n# runs:         " << runs;
  out << "\n# n_bins:       " << n_bins;
  out << "\n# temperature:  " << boost::any_cast<double>(settings["temperature"]);
  out << "\n# flatness:     " << boost::any_cast<double>(settings["flatness"]);
  out << "\n# one-over-t-c: " << boost::any_cast<double>(settings["one-over-t-c"]);
  out << "\n#" << std::endl;
}

/**
 * Public Methods
 */

/**
 * Constructor
 */
LennardJonesSystem::LennardJonesSystem(settings_t &s)
  : SimulationSystem(s), u_sphere_dist(3) {}

po::options_description LennardJonesSystem::get_program_options() {
  po::options_description desc("LennardJones System Options");
  desc.add_options()
    (
      "lj-size",
      po::value<double>(&size)->default_value(5),
      "Box length"
    ) (
      "lj-num-particles",
      po::value<size_t>(&num_particles)->default_value(2)->notifier(boost::bind(&LennardJonesSystem::set_particles, this, _1)),
      "Number of particles"
    ) (
      "lj-bins",
      po::value<size_t>(&n_bins)->default_value(100)->notifier(boost::bind(&LennardJonesSystem::set_bins, this, _1)),
      "Number of particles"
    ) (
      "lj-sigma",
      po::value<double>(&sigma)->default_value(1),
      "LJ sigma"
    ) (
      "lj-epsilon",
      po::value<double>(&epsilon)->default_value(1),
      "LJ epsilon"
    ) (
      "lj-cutoff-radius",
      po::value<double>(&cutoff_radius)->default_value(-1),
      "LJ cutoff radius, default: -1 (half of box size)"
    ) (
      "lj-emin",
      po::value<double>(&e_min)->default_value(-1),
      "LJ E_min"
    ) (
      "lj-emax",
      po::value<double>(&e_max)->default_value(3),
      "LJ E_max"
    ) (
      "lj-delta-r",
      po::value<double>(&delta_r)->default_value(1),
      "LJ Delta r"
    );
  return desc;
}

void LennardJonesSystem::parse_arguments(boost::program_options::variables_map &) {}

void LennardJonesSystem::setup() {
  SimulationSystem::setup();

  try {
    sampler       = boost::any_cast<size_t>(settings["sampler"]);
  } catch(const boost::bad_any_cast &e) {
    std::cerr << "Error while extracting settings in "
      << __FILE__ << ":" << __LINE__ << "\nExecption: "
      << e.what() << std::endl;
    throw e;
  }

  // set_size should have been called at this point
  // and thus n_bins should have a value
  settings["macrostates"] = n_bins;
  size_t num_unvisited_energies = 1;
  settings["num_unvisited_energies"] = num_unvisited_energies;

  //read_exact_dos();
  setup_output();

  // precalculate constants
  sigma6           = sigma*sigma*sigma*sigma*sigma*sigma;
  four_epsilon     = 4*epsilon;
  cutoff_radius_sq = cutoff_radius*cutoff_radius;

  box_dimensions[0] = box_dimensions[1] = box_dimensions[2] = size;
}

bool LennardJonesSystem::run() {
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

  write_output();

  out << "# last Q/Qd matrix:\n# " << Q << "\n# " << Qd << std::endl;

  return true;
}

/* vim: set ts=2 sw=2 sts=2 tw=0 expandtab :*/
