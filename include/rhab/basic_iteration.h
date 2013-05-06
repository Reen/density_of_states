#ifndef BASIC_ITERATION_H

#define BASIC_ITERATION_H

namespace rhab {
	template <class T>
	class basic_iteration {
	public:
		basic_iteration(std::size_t max_iter, double relative_eps)
			: i(0), max_iter_(max_iter), relative_eps_(relative_eps) {}

		inline bool converged(T iter1, T iter1end, T iter2, double & dist) {
			bool converged = true;

			// break after maximum number of iterations
			if (i == max_iter_) return true;

			for (;iter1 < iter1end;
					++iter1, ++iter2) {
				if (*iter1 > 0) {
					if (converged && fabs((*iter2)/(*iter1)-1.0) > relative_eps_) {
						dist = (*iter2)/(*iter1);
						converged = false;
						break;
					}
				}
			}
			return converged;
		}

		inline basic_iteration<T>& operator++() {
			++i;
			return (*this);
		}
		inline bool first() { return i == 0; }
		inline const std::size_t& iterations() const { return i; }
		inline double relative_tolerance() { return relative_eps_; }
		inline std::size_t max_iterations() { return max_iter_; }
	protected:
		std::size_t i;
		std::size_t max_iter_;
		double relative_eps_;
	};
}

#endif /* end of include guard: BASIC_ITERATION_H */
