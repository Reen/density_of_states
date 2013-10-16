#ifndef MISC_H_DIGGKQCX
#define MISC_H_DIGGKQCX

#include <boost/math/common_factor.hpp>

namespace rhab {

	template<class T>
	int get_gcd(const T &a) {
		int curgcd = a[0];
		for (size_t i = 1; i < a.size()-1; i++) {
			curgcd = boost::math::gcd(curgcd, a[i]);
			if (curgcd == 1) {
				break;
			}
		}
		return curgcd;
	}

}
#endif /* end of include guard: MISC_H_DIGGKQCX */

