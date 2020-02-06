/*
 * silo.hpp
 *
 *  Created on: Feb 6, 2020
 *      Author: dmarce1
 */

#ifndef OCTOTIGER_SILO_HPP_
#define OCTOTIGER_SILO_HPP_



#include <octotiger/full_state.hpp>

#include <hpx/include/plain_actions.hpp>
#include <hpx/include/serialization.hpp>

#include <silo.h>


struct silo_zone {
	full_state state;
	std::array<general_vect<fixed_real,NDIM>,NCHILD> nodes;
	template<class Arc>
	void serialize(Arc&& a, unsigned) {
		a & state;
		a & nodes;
	}
};


void silo_reset();
void silo_add_zones(const std::vector<silo_zone>&);

HPX_DECLARE_PLAIN_ACTION(silo_reset);
HPX_DECLARE_PLAIN_ACTION(silo_add_zones);


#endif /* OCTOTIGER_SILO_HPP_ */
