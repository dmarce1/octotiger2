/*
 * tree.cpp
 *
 *  Created on: Feb 5, 2020
 *      Author: dmarce1
 */

#include <octotiger/problems.hpp>
#include <octotiger/refinements.hpp>
#include <octotiger/silo.hpp>
#include <octotiger/tree.hpp>

HPX_REGISTER_COMPONENT(hpx::components::component<tree>, tree);

int tree::inx;
bool tree::global_time;
int tree::max_level;
real tree::cfl;
const int tree::bw = 1;
hpx::lcos::local::mutex tree::mtx;
std::vector<std::shared_ptr<super_array<full_state>>> tree::data_arrays_;
std::array<std::vector<std::shared_ptr<super_array<conserved>>>, NDIM> tree::flux_arrays_;

void tree::static_init() {
	const auto opts = options::get();
	inx = opts.grid_size;
	global_time = opts.global_time;
	cfl = opts.cfl;
	max_level = opts.max_level;
}

void tree::initialize() {
	for (int dim = 0; dim < NDIM; dim++) {
		const fixed_real div = fixed_real(2) / fixed_real(inx << level_);
		space_volume_.begin(dim) = fixed_real(index_volume_.begin(dim)) * div - fixed_real(1.0);
		space_volume_.end(dim) = fixed_real(index_volume_.end(dim)) * div - fixed_real(1.0);
	}
	dx_ = real(space_volume_.end(0) - space_volume_.begin(0)) / inx;
	dt_ = 0.0;
	{
		std::lock_guard<hpx::lcos::local::mutex> lock(mtx);
		if (level_ >= data_arrays_.size()) {
			data_arrays_.resize(level_ + 1);
			for (int dim = 0; dim < NDIM; dim++) {
				flux_arrays_[dim].resize(level_ + 1);
			}
		}
		if (data_arrays_[level_] == nullptr) {
			data_arrays_[level_] = std::make_shared<super_array<full_state>>();
			for (int dim = 0; dim < NDIM; dim++) {
				flux_arrays_[dim][level_] = std::make_shared<super_array<conserved>>();
			}
		}
	}
	state_ptr_ = data_arrays_[level_];
	const auto vol = index_volume_.expand(bw);
	state_ptr_->add_volume(vol);
	for (int dim = 0; dim < NDIM; dim++) {
		auto vol = index_volume_;
		vol.end()[dim]++;
		flux_ptr_[dim] = flux_arrays_[dim][level_];
		flux_ptr_[dim]->add_volume(vol);
	}
	refinement_flag = 0;
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

node_attr tree::get_node_attributes() const {
	node_attr a;
	a.leaf = is_leaf();
	a.space_volume = space_volume_;
	a.locality = hpx::find_here();
	return a;
}

void tree::find_family(hpx::id_type parent, hpx::id_type self, std::vector<hpx::id_type> neighbors) {
	parent_ = parent;
	self_ = self;
	neighbors_ = std::move(neighbors);
	std::array<hpx::future<std::vector<hpx::id_type>>, NSIBLING> ncfuts;
	std::array<std::vector<hpx::id_type>, NSIBLING> nchildren;
	if (!is_leaf()) {
		for (int si = 0; si < NSIBLING; si++) {
			if (neighbors_[si] != hpx::invalid_id) {
				ncfuts[si] = hpx::async<get_children_action>(neighbors_[si]);
			} else {
				ncfuts[si] = hpx::make_ready_future(std::vector<hpx::id_type>());
			}
		}
		for (int si = 0; si < NSIBLING; si++) {
			nchildren[si] = ncfuts[si].get();
			if (nchildren[si].empty()) {
				nchildren[si] = std::vector<hpx::id_type>(NCHILD, hpx::invalid_id);
			}
		}
		std::array<hpx::future<void>, NCHILD> cfuts;
		for (int ci = 0; ci < NCHILD; ci++) {
			std::vector<hpx::id_type> cneighbors(NSIBLING);
			for (int dim = 0; dim < NDIM; dim++) {
				const auto ci0 = ci ^ (1 << dim);
				if (((ci >> dim) & 1) == 0) {
					cneighbors[2 * dim + 0] = nchildren[2 * dim + 0][ci0];
					cneighbors[2 * dim + 1] = children_[ci0];
				} else {
					cneighbors[2 * dim + 0] = children_[ci0];
					cneighbors[2 * dim + 1] = nchildren[2 * dim + 1][ci0];
				}
			}
			cfuts[ci] = hpx::async<find_family_action>(children_[ci], self, children_[ci], cneighbors);
		}
		hpx::wait_all(cfuts.begin(), cfuts.end());
	}
}

std::vector<hpx::id_type> tree::get_children() const {
	return children_;
}

void tree::create_children() {
	assert(children_.empty());
	std::array<hpx::future<hpx::id_type>, NCHILD> futs;
	children_.resize(NCHILD);
	for (int ci = 0; ci < NCHILD; ci++) {
		volume<int> this_vol;
		for (int dim = 0; dim < NDIM; dim++) {
			if (((ci >> dim) & 1) == 0) {
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

std::vector<real> tree::get_prolong_con() {
	std::vector<real> data;

	for (auto I = index_volume_.begin(); I != index_volume_.end(); index_volume_.inc_index(I)) {
		const auto &U = (*state_ptr_)[I].U;
		for (int f = 0; f < NF; f++) {
			data.push_back(U[f]);
		}
	}
	return data;
}

void tree::set_con(const std::vector<real> &data) {
	int i = 0;

	for (auto I = index_volume_.begin(); I != index_volume_.end(); index_volume_.inc_index(I)) {
		auto &U = (*state_ptr_)[I].U;
		for (int f = 0; f < NF; f++) {
			U[f] = data[i++];
		}
	}
}

void tree::set_initial_conditions() {
	if (is_leaf()) {
		static const auto f = get_init_func();
		for (auto I = index_volume_.begin(); I != index_volume_.end(); index_volume_.inc_index(I)) {
			(*state_ptr_)[I].U = f(X(I));
			(*state_ptr_)[I].t = 0;
			(*state_ptr_)[I].dt = 0;
		}
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<set_initial_conditions_action>(children_[ci]);
		}
		hpx::wait_all(futs.begin(), futs.end());
	}
}

void tree::send_silo() {
	static const auto root_loc = hpx::find_all_localities()[0];
	if (is_leaf()) {
		std::vector<silo_zone> zones;

		for (auto I = index_volume_.begin(); I != index_volume_.end(); index_volume_.inc_index(I)) {
			silo_zone z;
			z.state = (*state_ptr_)[I];
			for (int ci = 0; ci < NCHILD; ci++) {
				for (int dim = 0; dim < NDIM; dim++) {
					if (((ci >> dim) & 1) == 0) {
						z.nodes[ci][dim] = X(I, dim) - fixed_real(0.5) * dx_;
					} else {
						z.nodes[ci][dim] = X(I, dim) + fixed_real(0.5) * dx_;
					}
				}
			}
			zones.push_back(z);
		}
		silo_add_zones_action()(root_loc, std::move(zones));
	} else {
		std::array<hpx::future<void>, NCHILD> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs[ci] = hpx::async<send_silo_action>(children_[ci]);
		}
		hpx::wait_all(futs.begin(), futs.end());
	}
}

bool tree::check_for_refine(fixed_real t) {
	bool rc = false;
	static const auto rfunc = get_refinement_function();
	if (is_leaf()) {
		if (max_level > level_) {
			for (auto I = index_volume_.begin(); I != index_volume_.end(); index_volume_.inc_index(I)) {
				gradient dW;
				for (int dim = 0; dim < NDIM; dim++) {
					auto Ip = I;
					auto Im = I;
					Ip[dim]++;
					Im[dim]--;
					dW[dim] = ((*state_ptr_)[Ip].W - (*state_ptr_)[Im].W) * 0.5;
				}
				assert(rfunc != nullptr);
				if (rfunc((*state_ptr_)[I].W, dW)) {
					rc = true;
					break;
				}
			}
			if (rc) {
				create_children();
			}
		}
	} else {
		std::array<hpx::future<bool>, NCHILD> cfuts;
		for (int ci = 0; ci < NCHILD; ci++) {
			cfuts[ci] = hpx::async<check_for_refine_action>(children_[ci], t);
		}
		hpx::wait_all(cfuts.begin(), cfuts.end());
		for (int ci = 0; ci < NCHILD; ci++) {
			rc = rc || cfuts[ci].get();
		}
	}
	return rc;
}
