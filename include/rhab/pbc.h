#ifndef PBC_H_IDBQB80S
#define PBC_H_IDBQB80S

#include <cmath>

namespace rhab {

/**
 * PBC - Periodic Boundary Conditions algorithms.
 *
 * Defines some static functions to operate on point_type points.
 *
 * @author René Haber
 */
class PBC
{
public:
  /**
   * PBC position with return.
   */
  template<class point_type, class box_size_type>
  static point_type position(point_type point, const box_size_type& size) {
    for (std::size_t i = 0; i < point.size(); i++) {
      if (point[i] > size[i]) {
        point[i] -= size[i];
      } else if (point[i] < 0.0) {
        point[i] += size[i];
      }
    }
    return point;
  }

  /**
   * PBC position with inplace changes.
   */
  template<class point_type, class box_size_type>
  static void positionIP(point_type &point, const box_size_type& size) {
    for (std::size_t i = 0; i < point.size(); i++) {
      if (point[i] > size[i]) {
        point[i] -= size[i];
      } else if (point[i] < 0.0) {
        point[i] += size[i];
      }
    }
  }

  /**
   * PBC separation with return.
   *
   * Calculates the Distance between two points in the SimulationBox.
   */
  template<class point_type, class box_size_type>
  static typename point_type::value_type
  separation(const point_type& point1, const point_type& point2,
             const box_size_type& size) {
    typename point_type::value_type acc(0.0);
    for (std::size_t i = 0; i < point1.size(); i++) {
      typename point_type::value_type diff = point1[i] - point2[i];
      typename point_type::value_type t = size[i];
      typename point_type::value_type t05 = 0.5*t;
      if (diff > t05) {
        diff -= t;
      } else if (diff < -t05) {
        diff += t;
      }

      acc += diff*diff;
    }
    return std::sqrt(acc);
  }

  /**
   * PBC separation squared with return.
   *
   * Calculates the Distance between two points in the SimulationBox but does not sqrt the result.
   */
  template<class point_type, class box_size_type>
  static typename point_type::value_type
  separationSquared(const point_type& point1, const point_type& point2,
                    const box_size_type& size) {
    typename point_type::value_type acc(0.0);
    for (std::size_t i = 0; i < point1.size(); i++) {
      typename point_type::value_type diff = point1[i] - point2[i];
      typename point_type::value_type t = size[i];
      typename point_type::value_type t05 = 0.5*t;
      if (diff > t05) {
        diff -= t;
      } else if (diff < -t05) {
        diff += t;
      }

      acc += diff*diff;
    }
    return acc;
  }
};

}; // namespace rhab

#endif /* end of include guard: PBC_H_IDBQB80S */

/* vim: set ts=2 sw=2 sts=2 tw=0 expandtab :*/
