/*
 * refinements.hpp
 *
 *  Created on: Feb 7, 2020
 *      Author: dmarce1
 */

#ifndef OCTOTIGER_REFINEMENTS_HPP_
#define OCTOTIGER_REFINEMENTS_HPP_


#include <octotiger/full_state.hpp>

#include <functional>

using refinement_func = std::function<bool(const primitive&, const gradient&)>;

refinement_func  get_refinement_function();


#endif /* OCTOTIGER_REFINEMENTS_HPP_ */
