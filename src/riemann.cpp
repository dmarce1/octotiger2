/*
 * riemann.cpp
 *
 *  Created on: Feb 5, 2020
 *      Author: dmarce1
 */

#include <octotiger/riemann.hpp>

#include <cassert>

conserved exact_riemann(const primitive &WL, const primitive &WR, int dim) {
	real dP;
	real P0;
	real err;
	real s0;
	conserved F;

	static const auto opts = options::get();
	const real fgamma = opts.fgamma;
	const auto gam1 = fgamma - 1.0;
	const auto gam2 = fgamma + 1.0;
	const auto gam3 = gam1 / (2.0 * fgamma);
	const auto gam4 = gam2 / (2.0 * fgamma);

	const auto &rhoL = WL.rho;
	const auto &rhoR = WR.rho;
	const auto &PR = WR.p;
	const auto &PL = WL.p;

	const auto UR = WR.to_con();
	const auto UL = WL.to_con();

	const auto AL = 2.0 / (gam2 * rhoL);
	const auto AR = 2.0 / (gam2 * rhoR);
	const auto BL = gam1 / gam2 * PL;
	const auto BR = gam1 / gam2 * PR;
	const auto &uL = WL.v[dim];
	const auto &uR = WR.v[dim];
	const auto aR = sqrt(fgamma * PR / rhoR);
	const auto aL = sqrt(fgamma * PL / rhoL);

	const auto sOL = uL + 2.0 * aL / gam1;
	const auto sOR = uR - 2.0 * aR / gam1;

	const auto fL = [&](real P) {
		if (P > PL) {
			return (P - PL) * sqrt(AL / (P + BL));
		} else {
			return 2.0 * aL / gam1 * (pow(P / PL, gam3) - 1);
		}
	};

	const auto fR = [&](real P) {
		if (P > PR) {
			return (P - PR) * sqrt(AR / (P + BR));
		} else {
			return 2.0 * aR / gam1 * (pow(P / PR, gam3) - 1);
		}
	};

	const auto dfLdP = [&](real P) {
		if (P > PL) {
			return sqrt(AL / (BL + P)) * (1.0 - (P - PL) / (2.0 * (BL + P)));
		} else {
			return 1.0 / (rhoL * aL) * pow(PL / P, gam4);
		}
	};

	const auto dfRdP = [&](real P) {
		if (P > PR) {
			return sqrt(AR / (BR + P)) * (1.0 - (P - PR) / (2.0 * (BR + P)));
		} else {
			return 1.0 / (rhoR * aR) * pow(PR / P, gam4);
		}
	};

	const auto f = [&](real P) {
		return fL(P) + fR(P) + uR - uL;
	};

	const auto dfdP = [&](real P) {
		return dfLdP(P) + dfRdP(P);
	};

	const auto WrL = [=]() {
		primitive W;
		const auto tmp = 2.0 / gam2 + (gam1 / gam2 / aL * uL);
		W.rho = rhoL * pow(tmp, 2.0 / gam1);
		W.v = WL.v;
		W.v[dim] = 2.0 / gam2 * (aL + gam1 / 2.0 * uL);
		W.p = PL * pow(2.0 / gam2 + (gam1 / gam2 / aL * uL), 2.0 * fgamma / gam1);
		return W;
	};

	const auto WrR = [=]() {
		primitive W;
		const auto tmp = 2.0 / gam2 - (gam1 / gam2 / aR * uR);
		W.rho = rhoR * pow(tmp, 2.0 / gam1);
		W.v = WR.v;
		W.v[dim] = 2.0 / gam2 * (-aR + gam1 / 2.0 * uR);
		W.p = PR * pow(2.0 / gam2 - (gam1 / gam2 / aR * uR), 2.0 * fgamma / gam1);
		return W;
	};

	const auto WO = [=]() {
		primitive W;
		W.p = W.rho = 0.0;
		W.v = vect(0.0);
		return W;
	};

	const auto W0L = [&]() {
		primitive W;
		real rho0L;
		if (P0 > PL) {
			const auto num = P0 / PL + gam1 / gam2;
			const auto den = (gam1 / gam2) * P0 / PL + 1.0;
			rho0L = rhoL * num / den;
		} else {
			rho0L = rhoL * pow(P0 / PL, 1.0 / fgamma);
		}
		W.p = P0;
		W.rho = rho0L;
		W.v = WL.v;
		W.v[dim] = s0;
		return W;
	};

	const auto W0R = [&]() {
		primitive W;
		real rho0R;
		if (P0 > PR) {
			const auto num = P0 / PR + gam1 / gam2;
			const auto den = (gam1 / gam2) * P0 / PR + 1.0;
			rho0R = rhoR * num / den;
		} else {
			rho0R = rhoR * pow(P0 / PR, 1.0 / fgamma);
		}
		W.p = P0;
		W.rho = rho0R;
		W.v = WR.v;
		W.v[dim] = s0;
		return W;
	};

	primitive Wi;
	if (sOL < sOR) {
		if (sOL > 0.0) {
			if (uL - aL > 0.0) {
				Wi = WL;
			} else if (sOL < 0.0) {
				Wi = WO();
			} else {
				Wi = WrL();
			}
		} else if (sOR < 0.0) {
			if (uR + aR < 0.0) {
				Wi = WR;
			} else if (sOR > 0.0) {
				Wi = WO();
			} else {
				Wi = WrR();
			}
		} else {
			Wi = WO();
		}
		F = Wi.to_flux(dim);
	} else {

		P0 = (PL + PR) / 2.0;
		if (P0 == 0.0) {
			P0 = 1.0e-6;
		}
		dP = (PR - PL) / 2.0;
		int iter = 0;
		real Plast;
		do {
			const auto g = f(P0);
			if (g == 0.0) {
				break;
			}
			const auto dgdp = dfdP(P0);
			assert(dgdp != 0.0);
			dP = -g / dgdp;
			real c0 = pow(0.1, double(iter) / double(1000));
			P0 = min(max(P0 + c0 * dP, P0 / 2.0), P0 * 2.0);
			if (iter > 1) {
				err = abs(P0 - Plast) / (P0 + Plast) / 2.0;
			} else {
				err = real::max();
			}
			iter++;
//				printf("%i %e %e %e %e %e %e %e %e\n", iter, c0, PL, PR, uL, uR, P0, dP, err);
			if (iter >= 1000) {
				printf("Riemann solver failed to converge\n");
				abort();
				;
			}
			Plast = P0;
		} while (err > 1.0e-12 && P0 > real::min());

		s0 = 0.5 * (uL + uR) + 0.5 * (fR(P0) - fL(P0));
		const auto QL = sqrt((P0 + BL) / AL);
		const auto QR = sqrt((P0 + BR) / AR);
		real sL, sR, sHL, sTL, sHR, sTR;

		if (s0 > 0.0) {
			if (P0 > PL) {
				sL = uL - QL / rhoL;
				if (sL > 0.0) {
					Wi = WL;
				} else {
					Wi = W0L();
				}
			} else {
				const auto a0L = aL * pow(P0 / PL, (fgamma - 1) / (2 * fgamma));
				sHL = uL - aL;
				sTL = s0 - a0L;
				assert(sHL <= sTL);
				if (sHL > 0.0) {
					Wi = WL;
				} else if (sTL < 0.0) {
					Wi = W0L();
				} else {
					Wi = WrL();
				}
			}
			F = Wi.to_flux(dim);
		} else {
			if (P0 > PR) {
				sR = uR + QR / rhoR;
				if (sR < 0.0) {
					Wi = WR;
				} else {
					Wi = W0R();
				}
			} else {
				const auto a0R = aR * pow(P0 / PR, (fgamma - 1) / (2 * fgamma));
				sHR = uR + aR;
				sTR = s0 + a0R;
	//			assert(sHR >= sTR);
				if (sHR < 0.0) {
					Wi = WR;
				} else if (sTR > 0.0) {
					Wi = W0R();
				} else {
					Wi = WrR();
				}
			}
			F = Wi.to_flux(dim);
		}
	}
	return F;
}

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
	FR.S[dim] += WR.p;
	FL.S[dim] += WL.p;
	FR.E += WR.p * WR.v[dim];
	FL.E += WL.p * WL.v[dim];
	F = (FR + FL) * 0.5;
	return F;

}

conserved riemann_solver(const primitive &WL, const primitive &WR, int dim) {
	return kurganov_tadmor(WL, WR, dim);
}
