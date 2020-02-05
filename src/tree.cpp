/*
 * tree.cpp
 *
 *  Created on: Feb 5, 2020
 *      Author: dmarce1
 */

#include <octotiger/math.hpp>
#include <octotiger/options.hpp>
#include <octotiger/riemann.hpp>
#include <octotiger/tree.hpp>

HPX_REGISTER_COMPONENT(hpx::components::component<tree>, tree);

int tree::inx;
bool tree::global_time;
fixed_real tree::cfl;
const int tree::bw = 1;
std::vector<std::shared_ptr<super_array<full_state>>> tree::data_arrays_;

void tree::static_init() {
	const auto opts = options::get();
	inx = opts.grid_size;
	global_time = opts.global_time;
	cfl = opts.cfl;
}

void tree::initialize() {
	for (int dim = 0; dim < NDIM; dim++) {
		const auto div = 1.0 / (inx << level_);
		space_volume_.begin(dim) = index_volume_.begin(dim) * div;
		space_volume_.end(dim) = index_volume_.end(dim) * div;
	}
	dx_ = real(space_volume_.end(0) - space_volume_.begin(0)) / inx;
	dt_ = 0.0;
	if (level_ >= data_arrays_.size()) {
		data_arrays_.resize(level_ + 1);
	}
	if (data_arrays_[level_] == nullptr) {
		data_arrays_[level_] = std::make_shared<super_array<full_state>>();
	}
	state_ptr_ = data_arrays_[level_];
	const auto vol = index_volume_.expand(bw);
	state_ptr_->add_volume(vol);
}

void tree::set_as_root() {
	for (int dim = 0; dim < NDIM; dim++) {
		index_volume_.begin(dim) = 0;
		index_volume_.end(dim) = inx;
	}
	t_ = 0.0;
	level_ = 0;
	initialize();
	for (auto I = index_volume_.begin(); I != index_volume_.end(); index_volume_.inc_index(I)) {
		(*state_ptr_)[I].t = (*state_ptr_)[I].dt = fixed_real(0.0);
	}

}

void tree::create_children() {
	assert(children_.empty());
	std::array<hpx::future<hpx::id_type>, NCHILD> futs;
	children_.resize(NCHILD);
	for (int ci = 0; ci < NCHILD; ci++) {
		volume<int> this_vol;
		for (int dim = 0; dim < NDIM; dim++) {
			if ((ci >> dim) & 1 == 0) {
				this_vol.begin(dim) = 2 * index_volume_.begin(dim);
				this_vol.end(dim) = index_volume_.begin(dim) + index_volume_.end(dim);
			} else {
				this_vol.begin(dim) = index_volume_.begin(dim) + index_volume_.end(dim);
				this_vol.end(dim) = 2 * index_volume_.end(dim);
			}
		}
		futs[ci] = hpx::new_<tree>(hpx::find_here(), this_vol, level_ + 1, t_);
	}
	for (int ci = 0; ci < NCHILD; ci++) {
		children_[ci] = futs[ci].get();
	}
}

tree::tree(const volume<int> &vol, int lev, fixed_real t) :
		index_volume_(vol), level_(lev), t_(t) {
	initialize();
}

void tree::con_to_prim(fixed_real t, fixed_real dt) {
	if (is_leaf()) {
		if (global_time || t + dt == t_ + dt_) {
			for (auto I = index_volume_.begin(); I != index_volume_.end(); index_volume_.inc_index(I)) {
				primitive &W = (*state_ptr_)[I].W;
				const conserved &U = (*state_ptr_)[I].U;
				W = U.to_prim();
				t_ += dt_;
				(*state_ptr_)[I].t = t_;
			}
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<con_to_prim_action>(children_[ci], t, dt);
		}
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
	}
}

fixed_real tree::timestep(fixed_real t) {
	fixed_real dt = fixed_real::max();
	if (is_leaf()) {
		if (t == t_) {
			for (auto I = index_volume_.begin(); I != index_volume_.end(); index_volume_.inc_index(I)) {
				primitive W = (*state_ptr_)[I].W;
				const auto vsig = W.signal_speed();
				dt_ = fixed_real(dx_ / vsig) * cfl;
				dt_ = dt_.nearest_log2();
				dt_ = min(dt_, dt_.next_bin() - t);
				dt = min(dt_, dt);
				(*state_ptr_)[I].dt = dt_;

			}
		}
	} else {
		std::array<hpx::future<fixed_real>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<timestep_action>(children_[ci], t);
		}
		for (int ci = 0; ci < NCHILD; ci++) {
			dt = min(dt, futs[ci].get());
		}
	}
	return dt;
}

void tree::physical_bc_primitive() {
	if (is_leaf()) {
		for (int dim = 0; dim < NDIM; dim++) {
			if (space_volume_.begin(dim) == fixed_real(0.0)) {
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
				bc_vol.begin(dim) = 1 << level_;
				bc_vol.end(dim) = (1 << level_) + 1;
				for (auto I = bc_vol.begin(); I != bc_vol.end(); bc_vol.inc_index((I))) {
					auto Im = I;
					Im[dim]--;
					(*state_ptr_)[I].W = (*state_ptr_)[Im].W;
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

void tree::update_con(fixed_real t, fixed_real dt) {
	if (is_leaf()) {
		primitive WR;
		primitive WL;
		for (auto I = index_volume_.begin(); I != index_volume_.end(); index_volume_.inc_index(I)) {
			for (int dim = 0; dim < NDIM; dim++) {
				const auto &IR = I;
				auto IL = I;
				IL[dim]--;
				auto &R = (*state_ptr_)[IR];
				auto &L = (*state_ptr_)[IL];
				if (L.t + L.dt == t + dt || R.t + R.dt == t + dt || global_time) {
					const auto this_dt = real(global_time ? dt : min(L.dt, R.dt));
					WR = R.W - (R.dW[dim] - R.dWdt() * this_dt) * 0.5;
					WL = L.W + (L.dW[dim] + L.dWdt() * this_dt) * 0.5;
					const auto F = riemann_solver(WL, WR, dim);
					const auto dU = F * (this_dt / dx_);
					R.U = R.U + dU;
					L.U = L.U - dU;
				}
			}
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<update_con_action>(children_[ci], t, dt);
		}
	}

}
