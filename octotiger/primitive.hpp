/*
 * primitive.hpp
 *
 *  Created on: Feb 2, 2020
 *      Author: dmarce1
 */

#ifndef OCTOTIGER_PRIMITIVE_HPP_
#define OCTOTIGER_PRIMITIVE_HPP_

#include <octotiger/options.hpp>
#include <octotiger/vect.hpp>

#define NF (NDIM + 2)

class conserved;

struct primitive: public general_vect<real, NF> {
	real &rho;
	real &p;
	general_vect<real, NDIM> &v;
	primitive();
	primitive& operator=(const general_vect<real, NF> &other);
	primitive& operator=(const primitive &other);
	conserved to_con() const;
	real sound_speed() const;
	real signal_speed() const;
};


using gradient = general_vect<primitive,NDIM>;

#endif /* OCTOTIGER_PRIMITIVE_HPP_ */
