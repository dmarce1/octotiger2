#include <octotiger/math.hpp>
#include <octotiger/riemann.hpp>
#include <octotiger/tree.hpp>

bool tree::adjust_dt(fixed_real t) {
	bool rc = false;
	if (is_leaf()) {
		if (t_ == t) {
			std::array<hpx::future<fixed_real>, NSIBLING> dtfuts;
			for (int si = 0; si < NSIBLING; si++) {
				if (neighbors_[si] != hpx::invalid_id) {
					dtfuts[si] = hpx::async<get_dt_action>(neighbors_[si]);
				} else {
					dtfuts[si] = hpx::make_ready_future(fixed_real::max() / fixed_real(2));
				}
			}
			for (int si = 0; si < NSIBLING; si++) {
				const auto max_dt = fixed_real(2) * dtfuts[si].get();
				if (dt_ > max_dt) {
					{
						std::lock_guard<hpx::lcos::local::mutex> lock(mtx_);
						dt_ = max_dt;
					}
					rc = true;
				}
			}
			load_times();
		}
	} else {
		std::array<hpx::future<bool>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<adjust_dt_action>(children_[ci], t);
		}
		for (int ci = 0; ci < NCHILD; ci++) {
			rc = futs[ci].get() || rc;
		}
	}
	return rc;
}

void tree::con_to_prim(fixed_real t, fixed_real dt) {
	if (is_leaf()) {
		if (global_time || (t + dt == t_ + dt_)) {
			t_ = t + dt;
			for (auto I = index_volume_.begin(); I != index_volume_.end(); index_volume_.inc_index(I)) {
				primitive &W = (*W_ptr_)[I];
				const conserved &U = (*U_ptr_)[I];
				W = U.to_prim();
				(*t_ptr_)[I] = t + dt;
			}
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<con_to_prim_action>(children_[ci], t, dt);
		}
		hpx::wait_all(futs.begin(), futs.end());
	}
}

void tree::gradients(fixed_real t) {
	if (is_leaf()) {
		if (global_time || t == t_) {

			for (auto I = index_volume_.begin(); I != index_volume_.end(); index_volume_.inc_index(I)) {
				const primitive W = (*W_ptr_)[I];
				for (int dim = 0; dim < NDIM; dim++) {
					index_type Ip = I;
					index_type Im = I;
					Ip[dim]++;
					Im[dim]--;
					const primitive &Wp = (*W_ptr_)[Ip];
					const primitive &Wm = (*W_ptr_)[Im];
					primitive &dW = (*dW_ptr_)[I][dim];
					for (int f = 0; f < NF; f++) {
						dW[f] = minmod(Wp[f] - W[f], W[f] - Wm[f]);
					}
				}
			}
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<gradients_action>(children_[ci], t);
		}
		hpx::wait_all(futs.begin(), futs.end());
	}
}

void tree::load_times() {
	for (auto I = index_volume_.begin(); I != index_volume_.end(); index_volume_.inc_index(I)) {
		(*dt_ptr_)[I] = dt_;
	}
}

fixed_real tree::timestep(fixed_real t) {
	if (is_leaf()) {
		if (t == t_ || global_time) {
			dt_ = fixed_real::max();
			for (auto I = index_volume_.begin(); I != index_volume_.end(); index_volume_.inc_index(I)) {
				primitive W = (*W_ptr_)[I];
				const auto vsig = W.signal_speed();
				fixed_real dt = real(dx_) / vsig * cfl;
				dt_ = min(dt_, dt);
			}
			dt_ = dt_.nearest_log2();
			dt_ = min(dt_, t.next_bin() - t);
			load_times();
		}
	} else {
		std::array<hpx::future<fixed_real>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<timestep_action>(children_[ci], t);
		}
		dt_ = fixed_real::max();
		for (int ci = 0; ci < NCHILD; ci++) {
			dt_ = min(dt_, futs[ci].get());
		}
	}
	return dt_;
}

void tree::physical_bc_primitive() {
	if (is_leaf()) {

		for (int dim = 0; dim < NDIM; dim++) {
			if (space_volume_.begin(dim) == fixed_real(-1.0)) {
				volume<int> bc_vol = index_volume_;
				bc_vol.begin(dim) = -1;
				bc_vol.end(dim) = 0;
				for (auto I = bc_vol.begin(); I != bc_vol.end(); bc_vol.inc_index((I))) {
					auto Ip = I;
					Ip[dim]++;
					auto &W = (*W_ptr_)[I];
					W = (*W_ptr_)[Ip];
					W.v[dim] = min(W.v[dim], real(0.0));
				}
			}
			if (space_volume_.end(dim) == fixed_real(1.0)) {
				volume<int> bc_vol = index_volume_;
				bc_vol.begin(dim) = inx << level_;
				bc_vol.end(dim) = (inx << level_) + 1;
				for (auto I = bc_vol.begin(); I != bc_vol.end(); bc_vol.inc_index((I))) {
					auto Im = I;
					Im[dim]--;
					auto &W = (*W_ptr_)[I];
					W = (*W_ptr_)[Im];
					W.v[dim] = max(W.v[dim], real(0.0));
				}

			}
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<physical_bc_primitive_action>(children_[ci]);
		}
		hpx::wait_all(futs.begin(), futs.end());
	}
}

void tree::physical_bc_gradient() {
	if (is_leaf()) {

		for (int dim = 0; dim < NDIM; dim++) {
			if (space_volume_.begin(dim) == fixed_real(-1.0)) {
				volume<int> bc_vol = index_volume_;
				bc_vol.begin(dim) = -1;
				bc_vol.end(dim) = 0;
				for (auto I = bc_vol.begin(); I != bc_vol.end(); bc_vol.inc_index((I))) {
					auto Ip = I;
					Ip[dim]++;
					for (int this_dim = 0; this_dim < NDIM; this_dim++) {
						(*dW_ptr_)[I][this_dim] = this_dim == dim ? general_vect<real, NF>(0) : (*dW_ptr_)[Ip][this_dim];
					}
					(*t_ptr_)[I] = (*t_ptr_)[Ip];
					(*dt_ptr_)[I] = (*dt_ptr_)[Ip];
				}
			}
			if (space_volume_.end(dim) == fixed_real(1.0)) {
				volume<int> bc_vol = index_volume_;
				bc_vol.begin(dim) = inx << level_;
				bc_vol.end(dim) = (inx << level_) + 1;
				for (auto I = bc_vol.begin(); I != bc_vol.end(); bc_vol.inc_index((I))) {
					auto Im = I;
					Im[dim]--;
					for (int this_dim = 0; this_dim < NDIM; this_dim++) {
						(*dW_ptr_)[I][this_dim] = this_dim == dim ? general_vect<real, NF>(0) : (*dW_ptr_)[Im][this_dim];
					}
					(*t_ptr_)[I] = (*t_ptr_)[Im];
					(*dt_ptr_)[I] = (*dt_ptr_)[Im];
				}

			}
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<physical_bc_gradient_action>(children_[ci]);
		}
		hpx::wait_all(futs.begin(), futs.end());
	}
}

void tree::update_con(fixed_real t, fixed_real dt) {
	if (is_leaf()) {
		for (int dim = 0; dim < NDIM; dim++) {
			auto flux_volume = index_volume_;
			flux_volume.end(dim)++;for
(			auto I = flux_volume.begin(); I != flux_volume.end(); flux_volume.inc_index(I)) {
				const auto &IR = I;
				auto IL = I;
				IL[dim]--;
				auto& LU = (*U_ptr_)[IL];
				auto& RU = (*U_ptr_)[IR];
				const auto& Lt = (*t_ptr_)[IL];
				const auto& Ldt = (*dt_ptr_)[IL];
				const auto& Rt = (*t_ptr_)[IR];
				const auto& Rdt = (*dt_ptr_)[IR];
				if ((Lt + Ldt == t + dt) || (Rt + Rdt == t + dt) || global_time) {
					const auto this_dt = global_time ? dt : min(Ldt, Rdt);
					const auto& F = (*flux_ptr_[dim])[I];
					const auto dU = F * (real(this_dt) / real(dx_));
					if (index_volume_.contains(IR)) {
						RU = RU + dU;
					}
					if (index_volume_.contains(IL)) {
						LU = LU - dU;
					}
				}
			}
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<update_con_action>(children_[ci], t, dt);
		}
		hpx::wait_all(futs.begin(),futs.end());
	}
}

void tree::compute_fluxes(fixed_real t, fixed_real dt) {
	if (is_leaf()) {
		primitive WR;
		primitive WL;
		for (int dim = 0; dim < NDIM; dim++) {
			auto flux_volume = index_volume_;
			if (neighbors_[2 * dim + 1] == hpx::invalid_id) {
				flux_volume.end(dim)++;}
			for (auto I = flux_volume.begin(); I != flux_volume.end(); flux_volume.inc_index(I)) {
				const auto &IR = I;
				auto IL = I;
				IL[dim]--;
				const auto &dWR = (*dW_ptr_)[IR];
				const auto &dWL = (*dW_ptr_)[IL];
				const auto& Lt = (*t_ptr_)[IL];
				const auto& Ldt = (*dt_ptr_)[IL];
				const auto& Rt = (*t_ptr_)[IR];
				const auto& Rdt = (*dt_ptr_)[IR];
				if ((Lt + Ldt == t + dt) || (Rt + Rdt == t + dt) || global_time) {
					const auto this_dt = global_time ? dt : min(Ldt, Rdt);
					WR = (*W_ptr_)[IR];
					WL = (*W_ptr_)[IL];
					WR = WR - dWR[dim] * 0.5;
					WL = WL + dWL[dim] * 0.5;
					WR = WR + WR.dWdt(dWR[dim]) * real(this_dt + t - Rt) * 0.5;
					WL = WL + WL.dWdt(dWL[dim]) * real(this_dt + t - Lt) * 0.5;
					auto &F = (*flux_ptr_[dim])[I];
					F = riemann_solver(WL, WR, dim);
				}
			}
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<compute_fluxes_action>(children_[ci], t, dt);
		}
		hpx::wait_all(futs.begin(), futs.end());
	}
}
