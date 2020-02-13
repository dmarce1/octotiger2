#include <octotiger/silo.hpp>
#include <octotiger/tree.hpp>

#include <hpx/hpx_init.hpp>

#include <fenv.h>

int hpx_main(int argc, char *argv[]) {

	options opts;

	opts.process_options(argc, argv);

	tree::static_init();

	auto root = hpx::new_<tree>(hpx::find_here()).get();
	tree::set_as_root_action()(root);
	tree::set_initial_conditions_action()(root);
	while (tree::check_for_refine_action()(root, 0.0)) {
		tree::set_initial_conditions_action()(root);
	}
	tree::find_family_action()(root, hpx::invalid_id, root, std::vector<hpx::id_type>(NSIBLING, hpx::invalid_id));
	tree::con_to_prim_action()(root, 0.0, 0.0);
	tree::physical_bc_primitive_action()(root);
	tree::gradients_action()(root, 0.0);
	silo_begin();
	tree::send_silo_action()(root);
	silo_end("X.0.silo", 0.0);

	fixed_real t = 0.0;
	int step = 0;
	int oiter = 0;
	printf("Starting main execution loop...\n");
	while (t < fixed_real(opts.tmax)) {
		fixed_real dt = tree::timestep_action()(root, t);
		while (tree::adjust_dt_action()(root, t)) {
			NULL;
		}
		printf("%i %e %e\n", step, (double) t, (double) dt);
		tree::physical_bc_gradient_action()(root);
		tree::compute_fluxes_action()(root, t, dt);
		tree::update_con_action()(root, t, dt);
		tree::con_to_prim_action()(root, t, dt);
		t += dt;
		step++;
		tree::physical_bc_primitive_action()(root);
		tree::gradients_action()(root, t);
		const double t0 = t;
		const double dt0 = dt;
		const std::int64_t oi = t0 / opts.output_freq;
		const std::int64_t last_oi = (t0 - dt0) / opts.output_freq;
		if (oi != last_oi) {
			std::string name = "X." + std::to_string(++oiter) + ".silo";
			silo_begin();
			tree::send_silo_action()(root);
			silo_end(name, t);
		}
	}
	printf("%i %e\n", step, (double) t);
	printf("Finished!\n");
	return hpx::finalize();
}

int main(int argc, char *argv[]) {
	feenableexcept(FE_DIVBYZERO);
	feenableexcept(FE_INVALID);
	feenableexcept(FE_OVERFLOW);
	std::vector<std::string> cfg = { "hpx.commandline.allow_unknown=1" };
	hpx::init(argc, argv, cfg);
}
