/*
 * refinements.cpp
 *
 *  Created on: Feb 7, 2020
 *      Author: dmarce1
 */

#include <octotiger/refinements.hpp>

bool unigrid_refine(const primitive&, const gradient&) {
	return true;
}

bool density_refine(const primitive& W, const gradient& dW) {
	bool rc = false;
	for (int dim = 0; dim < NDIM; dim++) {
		rc = rc || ((dW[dim].rho / W.rho) > 0.1);
	}
	return rc;
}

refinement_func get_refinement_function() {
	static const auto opts = options::get();
	refinement_func f;
	if (opts.refinement == "den") {
		f = density_refine;
	} else if( opts.refinement == "uni") {
		f = unigrid_refine;
	} else {
		printf("Uknown refinement type %s\n", opts.refinement.c_str());
		abort();
	}
	return f;
}
