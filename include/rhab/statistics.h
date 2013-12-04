#ifndef STATISTICS_H_0PD5PSM1
#define STATISTICS_H_0PD5PSM1

// Boost Accumulator
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

namespace rhab {

struct StepStatistics {
	size_t step;
	boost::accumulators::accumulator_set<
		double,
		boost::accumulators::stats<
			boost::accumulators::tag::min,
			boost::accumulators::tag::max,
			boost::accumulators::tag::mean,
			boost::accumulators::tag::variance
				> > err1;
	boost::accumulators::accumulator_set<
		double,
		boost::accumulators::stats<
			boost::accumulators::tag::min,
			boost::accumulators::tag::max,
			boost::accumulators::tag::mean,
			boost::accumulators::tag::variance
				> > err2;
	boost::accumulators::accumulator_set<
		double,
		boost::accumulators::stats<
			boost::accumulators::tag::min,
			boost::accumulators::tag::max,
			boost::accumulators::tag::mean,
			boost::accumulators::tag::variance
				> > err3;
	boost::accumulators::accumulator_set<
		double,
		boost::accumulators::stats<
			boost::accumulators::tag::min,
			boost::accumulators::tag::max,
			boost::accumulators::tag::mean,
			boost::accumulators::tag::variance
				> > err4;
};

}

#endif /* end of include guard: STATISTICS_H_0PD5PSM1 */

