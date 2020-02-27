/*
 * tree.hpp
 *
 *  Created on: Feb 5, 2020
 *      Author: dmarce1
 */

#ifndef OCTOTIGER_TREE_HPP_
#define OCTOTIGER_TREE_HPP_

#include <octotiger/fixed_real.hpp>
#include <octotiger/full_state.hpp>
#include <octotiger/super_array.hpp>

#include <hpx/include/components.hpp>

struct node_attr {
	bool leaf;
	volume<int> index_volume;
	hpx::id_type locality;
	template<class A>
	void serialize(A &arc, unsigned) {
		arc & leaf;
		arc & index_volume;
		arc & locality;
	}

};

class tree: public hpx::components::component_base<tree> {

	static int inx_;
	static bool global_time_;
	static real cfl_;
	static int max_level_;
	const static int bw_;

	static std::vector<std::shared_ptr<super_array<conserved>>> U_arrays_;
	static std::vector<std::shared_ptr<super_array<primitive>>> W_arrays_;
	static std::vector<std::shared_ptr<super_array<gradient>>> dW_arrays_;
	static std::vector<std::shared_ptr<super_array<fixed_real>>> t_arrays_;
	static std::vector<std::shared_ptr<super_array<fixed_real>>> dt_arrays_;
	static std::array<std::vector<std::shared_ptr<super_array<conserved>>>, NDIM> flux_arrays_;

	volume<int> index_volume_;
	volume<fixed_real> space_volume_;
	int level_;
	fixed_real dx_;
	fixed_real t_;
	fixed_real dt_;

	std::shared_ptr<super_array<primitive>> W_ptr_;
	std::shared_ptr<super_array<gradient>> dW_ptr_;
	std::shared_ptr<super_array<conserved>> U_ptr_;
	std::shared_ptr<super_array<fixed_real>> t_ptr_;
	std::shared_ptr<super_array<fixed_real>> dt_ptr_;
	std::array<std::shared_ptr<super_array<conserved>>, NDIM> flux_ptr_;

	hpx::id_type parent_;
	hpx::id_type self_;
	std::vector<hpx::id_type> children_;
	std::vector<node_attr> children_attr_;
	std::array<hpx::id_type, NSIBLING> neighbors_;
	std::array<node_attr, NSIBLING> neighbor_attr_;

	std::atomic<int> refinement_flag;

	mutable hpx::lcos::local::mutex mtx_;

	void initialize();

public:

	tree() = default;
	tree(const volume<int>&, int, fixed_real);

	inline bool is_leaf() const {
		return children_.empty();
	}

	inline fixed_real X(const index_type &I, int dim) {
		return (fixed_real(I[dim]) + fixed_real(0.5)) * dx_ - fixed_real(1.0);
	}

	inline vect X(const index_type &I) {
		vect x;
		for (int dim = 0; dim < NDIM; dim++) {
			x[dim] = double(X(I, dim));
		}
		return x;
	}

	inline bool is_physical(int face) const {
		const auto dim = face / 2;
		if (face % 2 == 0) {
			return space_volume_.begin(dim) == fixed_real(-1.0);
		} else {
			return space_volume_.end(dim) == fixed_real(+1.0);
		}
	}

	inline bool is_fine_amr(int face) {
		return (is_leaf() && neighbors_[face] == hpx::invalid_id && !is_physical(face));
	}

	inline bool is_coarse_amr(int face) {
		if (neighbors_[face] != hpx::invalid_id) {
			return (is_leaf() && !neighbor_attr_[face].leaf && !is_physical(face));
		} else {
			return false;
		}
	}

	void load_times();

	void set_as_root();
	HPX_DEFINE_COMPONENT_ACTION(tree, set_as_root);

	void create_children();

	static void static_init();

	bool adjust_dt(fixed_real);
	HPX_DEFINE_COMPONENT_ACTION(tree, adjust_dt);

	fixed_real get_dt() const {
		std::lock_guard<hpx::lcos::local::mutex> lock(mtx_);
		return dt_;
	}
	HPX_DEFINE_COMPONENT_ACTION(tree, get_dt);

	void set_initial_conditions();
	HPX_DEFINE_COMPONENT_ACTION(tree, set_initial_conditions);

	fixed_real timestep(fixed_real);
	HPX_DEFINE_COMPONENT_ACTION(tree,timestep);

	void con_to_prim(fixed_real, fixed_real);
	HPX_DEFINE_COMPONENT_ACTION(tree, con_to_prim);

	void gradients(fixed_real t);
	HPX_DEFINE_COMPONENT_ACTION(tree, gradients);

	void physical_bc_primitive();
	HPX_DEFINE_COMPONENT_ACTION(tree,physical_bc_primitive);

	void physical_bc_gradient();
	HPX_DEFINE_COMPONENT_ACTION(tree,physical_bc_gradient);

	void amr_bc_primitive();
	HPX_DEFINE_COMPONENT_ACTION(tree,amr_bc_primitive);

	void amr_bc_gradient();
	HPX_DEFINE_COMPONENT_ACTION(tree,amr_bc_gradient);

	void update_con(fixed_real, fixed_real);
	HPX_DEFINE_COMPONENT_ACTION(tree, update_con);

	void send_silo();
	HPX_DEFINE_COMPONENT_ACTION(tree, send_silo);

	node_attr get_node_attributes() const;
	HPX_DEFINE_COMPONENT_ACTION(tree, get_node_attributes);

	void find_family(hpx::id_type, hpx::id_type, std::vector<hpx::id_type>);
	HPX_DEFINE_COMPONENT_ACTION(tree,find_family);

	std::vector<hpx::id_type> get_children() const;
	HPX_DEFINE_COMPONENT_ACTION(tree, get_children);

	bool check_for_refine(fixed_real t);
	HPX_DEFINE_COMPONENT_ACTION(tree, check_for_refine);

	void compute_fluxes(fixed_real t, fixed_real dt);
	HPX_DEFINE_COMPONENT_ACTION(tree,compute_fluxes);

	std::vector<sub_array<primitive>> get_prim_from_niece(volume<int>) const;
	HPX_DEFINE_COMPONENT_ACTION(tree, get_prim_from_niece);

	sub_array<primitive> get_restricted_prim(const volume<int>&) const;
	HPX_DEFINE_COMPONENT_ACTION(tree, get_restricted_prim);

	sub_array<primitive> get_prim_from_aunt(volume<int>) const;
	HPX_DEFINE_COMPONENT_ACTION(tree, get_prim_from_aunt);

	sub_array<primitive> get_prolonged_prim(const volume<int>&) const;
	HPX_DEFINE_COMPONENT_ACTION(tree, get_prolonged_prim);

	std::vector<sub_array<gradient>> get_gradient_from_niece(volume<int>) const;
	HPX_DEFINE_COMPONENT_ACTION(tree, get_gradient_from_niece);

	sub_array<gradient> get_restricted_gradient(const volume<int>&) const;
	HPX_DEFINE_COMPONENT_ACTION(tree, get_restricted_gradient);

	sub_array<gradient> get_gradient_from_aunt(volume<int>) const;
	HPX_DEFINE_COMPONENT_ACTION(tree, get_gradient_from_aunt);

	sub_array<gradient> get_prolonged_gradient(const volume<int>&) const;
	HPX_DEFINE_COMPONENT_ACTION(tree, get_prolonged_gradient);

};

#endif /* OCTOTIGER_TREE_HPP_ */
