/*
 * problems.cpp
 *
 *  Created on: Feb 6, 2020
 *      Author: dmarce1
 */

#include <octotiger/options.hpp>
#include <octotiger/problems.hpp>

static auto sod(const vect &x) {
	conserved U;
	U.S = vect(0);
	if( x[0] > 0.5) {
		U.D = 1.0;
		U.E = 2.5;
	} else {
		U.D = 0.125;
		U.E = 0.25;
	}
	return U;
}

init_func get_init_func() {
	static const auto opts = options::get();
	if (opts.problem == "sod") {
		return sod;
	} else {
		printf("%s is not a known problem\n", opts.problem.c_str());
	}
}
