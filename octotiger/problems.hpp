/*
 * problems.hpp
 *
 *  Created on: Feb 6, 2020
 *      Author: dmarce1
 */

#ifndef OCTOTIGER_PROBLEMS_HPP_
#define OCTOTIGER_PROBLEMS_HPP_


#include <octotiger/conserved.hpp>

#include <functional>

using init_func = std::function<conserved(const vect&)>;

init_func get_init_func();


#endif /* OCTOTIGER_PROBLEMS_HPP_ */
