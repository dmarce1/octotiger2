/*
 * full_state.cpp
 *
 *  Created on: Feb 5, 2020
 *      Author: dmarce1
 */

#include <octotiger/full_state.hpp>

primitive full_state::dWdt() const {
	primitive dwdt;
	for (int f = 0; f < NF; f++) {
		dwdt[f] = 0.0;
	}
	for (int dim = 0; dim < NDIM; dim++) {
		for (int f = 0; f < NF; f++) {
			dwdt[f] -= W.v[dim] * dW[dim][f];
		}
		dwdt.rho -= W.rho * dW[dim].v[dim];
		dwdt.v[dim] -= dW[dim].p / W.rho;
		dwdt.p -= W.v[dim] * dW[dim].p;
	}
	return dwdt;
}

