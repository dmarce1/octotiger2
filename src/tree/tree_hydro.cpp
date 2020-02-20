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
				primitive &W = (*state_ptr_)[I].W;
				const conserved &U = (*state_ptr_)[I].U;
				W = U.to_prim();
				(*state_ptr_)[I].t = t + dt;
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
				const primitive W = (*state_ptr_)[I].W;
				for (int dim = 0; dim < NDIM; dim++) {
					index_type Ip = I;
					index_type Im = I;
					Ip[dim]++;
					Im[dim]--;
					const primitive &Wp = (*state_ptr_)[Ip].W;
					const primitive &Wm = (*state_ptr_)[Im].W;
					primitive &dW = (*state_ptr_)[I].dW[dim];
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
		(*state_ptr_)[I].dt = dt_;
	}
}

fixed_real tree::timestep(fixed_real t) {
	if (is_leaf()) {
		if (t == t_ || global_time) {
			dt_ = fixed_real::max();
			for (auto I = index_volume_.begin(); I != index_volume_.end(); index_volume_.inc_index(I)) {
				primitive W = (*state_ptr_)[I].W;
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
					auto &W = (*state_ptr_)[I].W;
					W = (*state_ptr_)[Ip].W;
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
					auto &W = (*state_ptr_)[I].W;
					W = (*state_ptr_)[Im].W;
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
					auto &S = (*state_ptr_)[I];
					S.t = (*state_ptr_)[Ip].t;
					S.dt = (*state_ptr_)[Ip].dt;
					S.dW = (*state_ptr_)[Ip].dW;
					S.dW[dim] = general_vect<real, NF>(0);
				}
			}
			if (space_volume_.end(dim) == fixed_real(1.0)) {
				volume<int> bc_vol = index_volume_;
				bc_vol.begin(dim) = inx << level_;
				bc_vol.end(dim) = (inx << level_) + 1;
				for (auto I = bc_vol.begin(); I != bc_vol.end(); bc_vol.inc_index((I))) {
					auto Im = I;
					Im[dim]--;
					auto &S = (*state_ptr_)[I];
					S.t = (*state_ptr_)[Im].t;
					S.dt = (*state_ptr_)[Im].dt;
					S.dW = (*state_ptr_)[Im].dW;
					S.dW[dim] = general_vect<real, NF>(0);
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
			flux_volume.end(dim)++;/**/

			for
(			auto I = flux_volume.begin(); I != flux_volume.end(); flux_volume.inc_index(I)) {
				const auto &IR = I;
				auto IL = I;
				IL[dim]--;
				auto &R = (*state_ptr_)[IR];
				auto &L = (*state_ptr_)[IL];
				if ((L.t + L.dt == t + dt) || (R.t + R.dt == t + dt) || global_time) {
					const auto this_dt = global_time ? dt : min(L.dt, R.dt);
					const auto& F = (*flux_ptr_[dim])[I];
					const auto dU = F * (real(this_dt) / real(dx_));
					if (index_volume_.contains(IR)) {
						R.U = R.U + dU;
					}
					if (index_volume_.contains(IL)) {
						L.U = L.U - dU;
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
				auto &R = (*state_ptr_)[IR];
				auto &L = (*state_ptr_)[IL];
				if ((L.t + L.dt == t + dt) || (R.t + R.dt == t + dt) || global_time) {
					const auto this_dt = global_time ? dt : min(L.dt, R.dt);
					WR = R.W;
					WL = L.W;
					WR = WR - R.dW[dim] * 0.5;
					WL = WL + L.dW[dim] * 0.5;
					WR = WR + R.dWdt() * real(this_dt + t - R.t) * 0.5;
					WL = WL + L.dWdt() * real(this_dt + t - L.t) * 0.5;
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
