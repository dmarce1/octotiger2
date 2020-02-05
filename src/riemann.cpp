/*
 * riemann.cpp
 *
 *  Created on: Feb 5, 2020
 *      Author: dmarce1
 */

#include <octotiger/riemann.hpp>

conserved kurganov_tadmor(const primitive &WL, const primitive &WR, int dim) {
	const conserved UL = WL.to_con();
	const conserved UR = WR.to_con();
	conserved FR;
	conserved FL;
	conserved F;
	const auto aR = WR.sound_speed() + abs(WR.v[dim]);
	const auto aL = WL.sound_speed() + abs(WL.v[dim]);
	const auto a = max(aR, aL);
	for (int f = 0; f < NF; f++) {
		FR[f] = UR[f] * (WR.v[dim] - a);
		FL[f] = UL[f] * (WL.v[dim] + a);
	}
	FR.P[dim] += WR.p;
	FL.P[dim] += WL.p;
	FR.E += WR.p * WR.v[dim];
	FL.E += WL.p * WL.v[dim];
	F = (FR + FL) * 0.5;
	return F;

}

conserved riemann_solver(const primitive &WL, const primitive &WR, int dim) {
	return kurganov_tadmor(WL, WR, dim);
}
