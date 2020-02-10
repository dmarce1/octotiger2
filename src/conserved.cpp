/*
 * conserved.cpp
 *
 *  Created on: Feb 5, 2020
 *      Author: dmarce1
 */

#include <octotiger/conserved.hpp>

conserved& conserved::operator=(const conserved &other) {
	general_vect<real, NF>::operator=(other);
	return *this;
}

conserved& conserved::operator=(const general_vect<real, NF> &other) {
	general_vect<real, NF>::operator=(other);
	return *this;
}

conserved::conserved() :
		D((*this)[0]), E((*this)[1]), S(*(reinterpret_cast<general_vect<real, NDIM>*>(&((*this)[2])))) {
}

primitive conserved::to_prim() const {
	static const auto fgamma = options::get().fgamma;
	primitive W;
	W.rho = D;
	W.p = real(fgamma - 1.0) * (E - real(0.5) * S.dot(S) / D);
	W.v = S / D;
	return W;
}
