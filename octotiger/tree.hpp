/*
 * tree.hpp
 *
 *  Created on: Feb 5, 2020
 *      Author: dmarce1
 */

#ifndef OCTOTIGER_TREE_HPP_
#define OCTOTIGER_TREE_HPP_

#include <octotiger/conserved.hpp>
#include <octotiger/fixed_real.hpp>
#include <octotiger/primitive.hpp>
#include <octotiger/super_array.hpp>

#include <hpx/include/components.hpp>

struct full_state {
	conserved U;
	primitive W;
	gradient dW;
	fixed_real t;
	fixed_real dt;
};

class tree: public hpx::components::component_base<tree> {

	static int inx;
	static bool global_time;
	static fixed_real cfl;
	const static int bw;
	static std::vector<std::shared_ptr<super_array<full_state>>> data_arrays_;

	volume<int> index_volume_;
	volume<fixed_real> space_volume_;
	int level_;
	real dx_;
	fixed_real t_;
	fixed_real dt_;

	std::shared_ptr<super_array<full_state>> state_ptr_;

	std::vector<hpx::id_type> children_;

	void initialize();

public:

	tree() = default;
	tree(const volume<int>&, int, fixed_real);

	inline bool is_leaf() const {
		return children_.empty();
	}

	void set_as_root();

	void create_children();

	static void static_init();

	fixed_real timestep(fixed_real);
	HPX_DEFINE_COMPONENT_ACTION(tree,timestep);

	void con_to_prim(fixed_real, fixed_real);
	HPX_DEFINE_COMPONENT_ACTION(tree, con_to_prim);

	void gradients(fixed_real t);
	HPX_DEFINE_COMPONENT_ACTION(tree, gradients);

	void physical_bc_primitive();
	HPX_DEFINE_COMPONENT_ACTION(tree,physical_bc_primitive);

	void update_con(fixed_real, fixed_real);
	HPX_DEFINE_COMPONENT_ACTION(tree, update_con);
};

#endif /* OCTOTIGER_TREE_HPP_ */
