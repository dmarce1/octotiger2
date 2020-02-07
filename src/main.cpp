
#include <octotiger/silo.hpp>
#include <octotiger/tree.hpp>

#include <hpx/hpx_init.hpp>

int hpx_main(int argc, char *argv[]) {
	options opts;

	opts.process_options(argc,argv);

	tree::static_init();

	auto root = hpx::new_<tree>(hpx::find_here()).get();
	tree::set_as_root_action()(root);
	tree::set_initial_conditions_action()(root);
	tree::con_to_prim_action()(root, 0.0, 0.0);
	tree::physical_bc_primitive_action()(root);
	tree::gradients_action()(root, 0.0);
	tree::physical_bc_gradient_action()(root);
	silo_begin();
	tree::send_silo_action()(root);
	silo_end( "X.0.silo", 0.0);


	return hpx::finalize();
}

int main(int argc, char *argv[]) {
	std::vector<std::string> cfg = { "hpx.commandline.allow_unknown=1" };
	hpx::init(argc, argv, cfg);
}
