/*
 * silo.cpp
 *
 *  Created on: Feb 6, 2020
 *      Author: dmarce1
 */
#include <octotiger/silo.hpp>
#include <octotiger/fixed_real.hpp>

#include <hpx/lcos/local/mutex.hpp>

#include <unordered_map>

using node_type = general_vect<fixed_real, NDIM>;

struct node_hash {
	std::size_t operator()(const node_type &n) const {
		std::size_t key = 0;
		for (int dim = 0; dim < NDIM; dim++) {
			key ^= n[dim].get_int();
		}
		return key;
	}
};

static std::vector<silo_zone> zones_;
static hpx::lcos::local::mutex mtx_;
static std::unordered_map<node_type, int, node_hash> node_map_;
static int next_node_num_;

HPX_PLAIN_ACTION(silo_reset);
HPX_PLAIN_ACTION(silo_add_zones);

void silo_reset() {
	zones_.clear();
	node_map_.clear();
	next_node_num_ = 0;
}

void silo_add_zones(const std::vector<silo_zone> &zones) {
	std::lock_guard<hpx::lcos::local::mutex> lock(mtx_);

	for (const auto &z : zones) {
		for (int ci = 0; ci < NCHILD; ci++) {
			const auto iter = node_map_.find(z.nodes[ci]);
			if (iter == node_map_.end()) {
				node_map_[z.nodes[ci]] = next_node_num++;
			}
		}
	}

}

