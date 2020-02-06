/*
 * full_state.hpp
 *
 *  Created on: Feb 5, 2020
 *      Author: dmarce1
 */

#ifndef OCTOTIGER_FULL_STATE_HPP_
#define OCTOTIGER_FULL_STATE_HPP_

#include <octotiger/conserved.hpp>
#include <octotiger/fixed_real.hpp>
#include <octotiger/primitive.hpp>

struct full_state {
	conserved U;
	primitive W;
	gradient dW;
	fixed_real t;
	fixed_real dt;
	primitive dWdt() const;
	template<class A>
	void serialize(A &a, unsigned) {
		a & U;
		a & W;
		a & dW;
		a & t;
		a & dt;
	}
};

#endif /* OCTOTIGER_FULL_STATE_HPP_ */
