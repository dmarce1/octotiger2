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
	if (x[0] > 0.0) {
		U.D = 1.0;
		U.E = 2.5;
	} else {
		U.D = 0.125;
		U.E = 0.25;
	}
	return U;
}

static auto blast(const vect &x) {
	conserved U;
	U.S = vect(0);
	U.D = 1.0;
	U.E = max(1.0e-10, exp(-500.0 * x.dot(x)));
	return U;
}

init_func get_init_func() {
	static const auto opts = options::get();
	if (opts.problem == "sod") {
		return sod;
	} else if (opts.problem == "blast") {
		return blast;
	} else {
		printf("%s is not a known problem\n", opts.problem.c_str());
	}
}
