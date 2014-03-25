
#ifndef ACCUMULATOR_H_S1IRKODK
#define ACCUMULATOR_H_S1IRKODK
#include <algorithm>
#include <cstddef>
#include <limits>

namespace rhab {

class Accumulator {
  size_t count_;
  double mean_acc_;
  double var_acc_;

public:
  Accumulator() : count_(0), mean_acc_(0.), var_acc_(0.) {}

  void operator()(double value) {
    count_++;
    mean_acc_+= value;
    var_acc_ += value*value;
  }

  void operator+=(const Accumulator& other) {
    count_    += other.count_;
    mean_acc_ += other.mean_acc_;
    var_acc_  += other.var_acc_;
  }

  const size_t & count() const {
    return count_;
  }

  double mean() const {
    return mean_acc_/count_;
  }
};

class AccumulatorMinMax {
  size_t count_;
  double mean_acc_;
  double var_acc_;
  double min_;
  double max_;

public:
  AccumulatorMinMax()
    : count_(0), mean_acc_(0.), var_acc_(0.),
      min_(std::numeric_limits<double>::max()),
      max_(std::numeric_limits<double>::min())
  {}

  void operator()(double value) {
    count_++;
    mean_acc_+= value;
    var_acc_ += value*value;
    if (value < min_) {
      min_ = value;
    }
    if (value > max_) {
      max_ = value;
    }
  }

  void operator+=(const AccumulatorMinMax& other) {
    count_    += other.count_;
    mean_acc_ += other.mean_acc_;
    var_acc_  += other.var_acc_;
    min_       = std::min(min_, other.min_);
    max_       = std::max(max_, other.max_);
  }

  const size_t & count() const {
    return count_;
  }

  // implement better online algorithm found here:
  // http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
  double variance() const {
    return (var_acc_ - (mean_acc_*mean_acc_)/count_)/(count_-1);
  }

  double mean() const {
    return mean_acc_/count_;
  }

  const double& max() const {
    return max_;
  }

  const double& min() const {
    return min_;
  }
};

}


#endif /* end of include guard: ACCUMULATOR_H_S1IRKODK */

/* vim: set ts=2 sw=2 sts=2 tw=0 expandtab :*/
