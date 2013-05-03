#ifndef ISING_SYSTEM_H

#define ISING_SYSTEM_H

// C++ Standard Library
#include <fstream>
#include <iomanip>
#include <vector>

// Boost Accumulator
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/variance.hpp>

// Boost Assert
#include <boost/assert.hpp>

// Boost Bind
#include <boost/bind.hpp>

// Boost Format
#include <boost/format.hpp>

// Boost Program Options
#include <boost/program_options.hpp>

// Boost Random
#include <boost/random/uniform_01.hpp>


#include "simulation_system.h"
#include "mc_sampler.h"
#include "q_matrix_tools.h"

namespace po = boost::program_options;

using namespace rhab;

struct StepStatistics {
	size_t step;
	boost::accumulators::accumulator_set<
		double,
		boost::accumulators::stats<
			boost::accumulators::tag::min,
			boost::accumulators::tag::max,
			boost::accumulators::tag::mean,
			boost::accumulators::tag::variance,
			boost::accumulators::tag::median
				> > err1;
	boost::accumulators::accumulator_set<
		double,
		boost::accumulators::stats<
			boost::accumulators::tag::min,
			boost::accumulators::tag::max,
			boost::accumulators::tag::mean,
			boost::accumulators::tag::variance,
			boost::accumulators::tag::median
				> > err2;
};

class IsingSystem : public SimulationSystem {
private:
	typedef boost::numeric::ublas::matrix<signed char> storage_t;
	typedef std::vector< StepStatistics > error_acc_t;

	// size L of one dimension of the lattice
	size_t size;

	// the lattice itself
	storage_t lattice;

	/**
	 * Settings
	 */
	// the sampler chosen
	size_t sampler;

	// the interaction constant J
	double J;


	/**
	 * Misc
	 */
	boost::uniform_01<> dist01;

	// error accumulator array
	error_acc_t error_acc;

	// exact dos
	vector_double_t dos_exact_norm;

	// size dependent constants:
	int e_min;
	int e_max;
	size_t n_bins;

	void set_size(size_t L) {
		size = L;
		lattice.resize(L, L);
		e_min = -2 * L * L;
		//e_max =  2 * L * L;
		e_max =  0;
		n_bins = (e_max - e_min)/4 + 1;
		Q.resize(n_bins, n_bins);
		Qd.resize(n_bins, n_bins);
	}

	int reset() {
		int magnetization = 0;
		// set lattice to random state
		for (size_t i = 0; i < lattice.data().size(); ++i) {
			lattice.data()[i] = (dist01(rng.rng) > 0.5 ? 1 : -1);
			magnetization += lattice.data()[i];
		}
		return magnetization;
	}

	int calculate_energy() {
		int energy = 0;
		for (size_t i = 0; i < size; ++i) {
			for (size_t j = 0; j < size; ++j) {
				energy += -J * (lattice(i, j) * (
							lattice((i-1+size)%size,  j              ) +
							lattice( i             , (j-1+size)%size )
							));
			}
		}
		return energy;
	}

	template<class Sampler>
	void mc_loop() {
		
		for (size_t run = 0; run < runs; run++) {
			// variables for error / statistics calculation
			size_t error_check_freq = error_check_f;
			size_t index = 0;
			// reset Q matrix
			Q *= 0;
			int magnetization(0), energy(1);
			while (energy > 0) {
				magnetization = reset();
				energy        = calculate_energy();
				(void)magnetization;
			}
			int i, j;

			Sampler sampler(rng.rng, Q, settings);
			for (size_t step = 1; step <= steps; step++) {
				// single spin flip:
				i = (int)(dist01(rng.rng) * size);
				j = (int)(dist01(rng.rng) * size);

				// calculate change in energy
				int dE = 2 * J * lattice(i,j) * (
						lattice((i+1)%size     ,  j             ) +
						lattice((i-1+size)%size,  j             ) +
						lattice( i             , (j+1)%size     ) +
						lattice( i             , (j-1+size)%size)
						);

				BOOST_ASSERT((energy-e_min)%4 == 0);
				BOOST_ASSERT((energy-e_min+dE)%4 == 0);
				int Ei = (energy-e_min) / 4;
				int Ej = (energy-e_min+dE) / 4;
				BOOST_ASSERT(Ei >= 0);
				BOOST_ASSERT(Ej >= 0);

				bool out_of_bounds = false;
				// update Q matrix
				if (Ej < 0 || Ej >= (int)n_bins) {
					Q(Ei, Ei)++;
					out_of_bounds = true;
					sampler.rejected(Ei);
				} else {
					Q(Ei, Ej)++;
				}

				if (!out_of_bounds && sampler(energy, energy+dE, Ei, Ej)) {
					lattice(i, j) *= -1;
					energy += dE;
				}
				if (step % error_check_freq == 0) {
					// calculate_statistics()
					error_acc[index].step = step;
					Qd = normalize_q(Q);

					double err = calculate_error_q(dos_exact_norm, Qd);
					error_acc[index].err1(err);
					if (sampler.has_own_statistics()) {
						err = sampler.calculate_error(dos_exact_norm);
						error_acc[index].err2(err);
					}

					index++;
					if (step % (10*error_check_freq) == 0) {
						error_check_freq *= 10;
					}
				}
				sampler.check(step, run);
			}
		}
	}

	void read_exact_dos() {
		std::string format = "%1%/data/ising/exact_%2%.dat";
		std::string fn = str( boost::format(format)
				% boost::any_cast<std::string>(settings["executable_path"])
				% size
				);
		std::ifstream in(fn.c_str());

		dos_exact_norm.resize(n_bins);
		dos_exact_norm *= 0;
		for (size_t i = 0; i < n_bins; i++) {
			int energy;
			double dos;
			in >> energy >> dos;
			dos_exact_norm[i] = dos;
		}
		//std::cout << e_min << " " << e_max << " " << n_bins << std::endl;
		//std::cout << dos_exact_norm << std::endl;
		//std::exit(EXIT_FAILURE);
	}

	std::ofstream out;
	std::string tag;
	void setup_output() {
		std::vector<std::string> sampler_string;
		// initialize list of available samplers
		sampler_string.push_back("BM");
		sampler_string.push_back("WL");
		sampler_string.push_back("QB");
		sampler_string.push_back("QA");

		std::string format;
		switch (sampler) {
			case 1:
				format = "dos_%1%_%2%S_%3%R_%4%M_%6$0.2ff%7%%8%.out";
				break;
			case 2:
			case 3:
				format = "dos_%1%_%2%S_%3%R_%4%M%7%%8%.out";
				break;
			case 0:
			default:
				format = "dos_%1%_%2%S_%3%R_%4%M_%5$0.2fT%7%%8%.out";
		}
		std::string buf = str( boost::format(format)
				% sampler_string[sampler]
				% steps
				% runs
				% n_bins
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

public:
	IsingSystem()
	{}

	virtual boost::program_options::options_description get_program_options() {
		boost::program_options::options_description desc("Ising System Options");
		desc.add_options()
			("ising-size",        po::value<size_t>(&size)->default_value(14)->notifier(boost::bind(&IsingSystem::set_size, this, _1)), "Number of Spins in one dimension.\nGrid will be size*size")
			("ising-interaction", po::value<double>(&J)->default_value(1),      "Interaction J\nJ>0 - ferromagnetic\nJ<0 - antiferromagnetic\nJ=0 - noninteracting")
			;
		return desc;
	}

	virtual void setup(settings_t s) {
		settings = s;

		set_size(size);

		if (settings.count("seed")) {
			rng.seed = boost::any_cast<size_t>(settings["seed"]);
			rng.seed_set = true;
		} else {
			rng.initialize();
		}

		try {
			runs          = boost::any_cast<size_t>(settings["runs"]);
			steps         = boost::any_cast<size_t>(settings["steps"]);
			error_check_f = boost::any_cast<size_t>(settings["error_check_f"]);
			sampler       = boost::any_cast<size_t>(settings["sampler"]);
			tag           = boost::any_cast<std::string>(settings["tag"]);
		} catch(const boost::bad_any_cast &e) {
			std::cerr << "bla: " << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}

		// safety check
		if (steps % error_check_f != 0) {
			std::cerr << "Error: check-frequency and number of steps don't match" << std::endl;
			std::exit(EXIT_FAILURE);
		}
		size_t error_acc_size = 9 * (size_t)log10(steps / error_check_f) + 1;
		error_acc.resize(error_acc_size);

		settings["macrostates"] = n_bins;

		read_exact_dos();
		setup_output();
	}

	virtual bool run() {
		switch (sampler) {
		case 0:
			mc_loop<BoltzmannSampler>();
			break;
		case 1:
			mc_loop<WangLandauSampler>();
			break;
		/*case 2:
			mc_loop<QualityMeasureBSampler>();
			break;
		case 3:
			mc_loop<QualityMeasureASampler>();
			break;*/
		default:
			std::cerr << "Error: unknown sampler" << std::endl;
			return false;
		}

		out << "#\n#          time     mean_error      var_error   mean_error_s    var_error_s       qm_error\n";
		for (size_t i = 0; i < error_acc.size(); i++) {
			size_t time  = error_acc[i].step;
			double mean  = boost::accumulators::mean(error_acc[i].err1);
			double var   = boost::accumulators::variance(error_acc[i].err1);
			double mean2 = boost::accumulators::mean(error_acc[i].err2);
			double var2  = boost::accumulators::variance(error_acc[i].err2);
			//double qm_err   = error_acc[i].get<5>() / runs;
			out << std::setw(15) << std::right
				<< time
				<< std::setw(15) << std::right
				<< mean
				<< std::setw(15) << std::right
				<< var
				<< std::setw(15) << std::right
				<< mean2
				<< std::setw(15) << std::right
				<< var2
				//<< std::setw(15) << std::right
				//<< qm_err
				<< "\n";

		}
		return true;
	}
};

#endif /* end of include guard: ISING_SYSTEM_H */
