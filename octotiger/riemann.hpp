/*
 * riemann.hpp
 *
 *  Created on: Feb 2, 2020
 *      Author: dmarce1
 */

#ifndef OCTOTIGER_RIEMANN_HPP_
#define OCTOTIGER_RIEMANN_HPP_

#include <octotiger/conserved.hpp>
#include <octotiger/primitive.hpp>


conserved riemann_solver(const primitive &WL, const primitive &WR, int dim);

#endif /* OCTOTIGER_RIEMANN_HPP_ */
