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


void silo_begin();
void silo_end(const std::string&, fixed_real);
void silo_add_zones(const std::vector<silo_zone>&);

HPX_DEFINE_PLAIN_ACTION(silo_add_zones);


#endif /* OCTOTIGER_SILO_HPP_ */
