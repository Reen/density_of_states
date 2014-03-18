
#ifndef ACCUMULATOR_H_S1IRKODK
#define ACCUMULATOR_H_S1IRKODK
#include <cstddef>

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

}


#endif /* end of include guard: ACCUMULATOR_H_S1IRKODK */

/* vim: set ts=2 sw=2 sts=2 tw=0 expandtab :*/
