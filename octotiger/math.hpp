/*
 * math.hpp
 *
 *  Created on: Dec 6, 2019
 *      Author: dmarce1
 */

#ifndef SRC_MATH_HPP_
#define SRC_MATH_HPP_

#include <octotiger/vect.hpp>

vect rotate_to(const vect &u, const vect &n);
vect rotate_from(const vect &u, vect n);
void enable_floating_point_exceptions();

template<class T>
inline T minmod(const T &a, const T &b) {
	return (copysign(T(0.5), a) + copysign(T(0.5), b)) * min(abs(a), abs(b));
}

#endif /* SRC_MATH_HPP_ */
