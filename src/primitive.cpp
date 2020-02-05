/*
 * primitive.cpp
 *
 *  Created on: Feb 5, 2020
 *      Author: dmarce1
 */

#include <octotiger/conserved.hpp>
#include <octotiger/primitive.hpp>



primitive& primitive::operator=(const primitive &other) {
	general_vect<real, NF>::operator=(other);
	return *this;
}

primitive::primitive() :
		rho((*this)[0]), p((*this)[1]), v(*(reinterpret_cast<general_vect<real, NDIM>*>(&((*this)[2])))) {
}

conserved primitive::to_con() const {
	static const auto fgamma = options::get().fgamma;
	conserved U;
	U.rho = rho;
	U.P = v * rho;
	U.E = (p / real(fgamma - 1.0)) + real(0.5) * v.dot(v) * rho;
	return U;
}

real primitive::sound_speed() const {
	static const auto fgamma = options::get().fgamma;
	return sqrt(max(p, real(0)) / rho);
}

real primitive::signal_speed() const {
	static const auto fgamma = options::get().fgamma;
	return sound_speed() + abs(v);
}
