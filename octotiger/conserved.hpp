/*
 * conserved.hpp
 *
 *  Created on: Feb 2, 2020
 *      Author: dmarce1
 */

#ifndef OCTOTIGER_CONSERVED_HPP_
#define OCTOTIGER_CONSERVED_HPP_

#include <octotiger/options.hpp>
#include <octotiger/primitive.hpp>


struct conserved: public general_vect<real, NF> {
	real &rho;
	real &E;
	general_vect<real, NDIM> &P;
	conserved();
	conserved& operator=(const conserved& other);
	conserved& operator=(const general_vect<real, NF>&);
	primitive to_prim() const;

};


#endif /* OCrealOTIGER_CONSERVED_HPP_ */
