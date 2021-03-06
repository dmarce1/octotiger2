#include <octotiger/options.hpp>

#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/async.hpp>
#include <hpx/runtime/threads/run_as_os_thread.hpp>

#include <boost/program_options.hpp>

#include <fstream>
#include <iostream>
#include <thread>

options options::global;

HPX_PLAIN_ACTION(options::set, set_options_action);

options& options::get() {
	return global;
}

void options::set(options o) {
	global = o;
}

bool options::process_options(int argc, char *argv[]) {
	namespace po = boost::program_options;

	po::options_description command_opts("options");

	command_opts.add_options() //
	("cfl", po::value<double>(&cfl)->default_value(0.4), "CFL factor") //
	("config_file", po::value<std::string>(&config_file)->default_value(""), "configuration file") //
	("fgamma", po::value<double>(&fgamma)->default_value(7.0 / 5.0), "gamma for fluid gamma law") //
	("global_time", po::value<bool>(&global_time)->default_value(true), "enable global time-stepping") //
	("grid_size", po::value<int>(&grid_size)->default_value(32), "size of grid") //
	("help", "produce help message") //
	("max_level", po::value<int>(&max_level)->default_value(3), "size of grid") //
	("problem", po::value<std::string>(&problem)->default_value("sod"), "problem - blast, sod") //
	("output_freq", po::value<double>(&output_freq)->default_value(1.0), "output frequency") //
	("refinement", po::value<std::string>(&refinement)->default_value("den"), "refinement type - den, uni") //
	("tmax", po::value<double>(&tmax)->default_value(1.0), "tmax") //
			;

	boost::program_options::variables_map vm;
	po::store(po::parse_command_line(argc, argv, command_opts), vm);
	po::notify(vm);
	if (vm.count("help")) {
		std::cout << command_opts << "\n";
		return false;
	}
	if (!config_file.empty()) {
		std::ifstream cfg_fs { vm["config_file"].as<std::string>() };
		if (cfg_fs) {
			po::store(po::parse_config_file(cfg_fs, command_opts), vm);
		} else {
			printf("Configuration file %s not found!\n", config_file.c_str());
			return false;
		}
	}
	po::notify(vm);

	set(*this);

#define SHOW( opt ) std::cout << std::string( #opt ) << " = " << std::to_string(opt) << '\n';
#define SHOW_STRING( opt ) std::cout << std::string( #opt ) << " = " << opt << '\n';
	SHOW(cfl);
	SHOW_STRING(config_file);
	SHOW(fgamma);
	SHOW(global_time);
	SHOW(grid_size);
	SHOW(max_level);
	SHOW(output_freq);
	SHOW_STRING(problem);
	SHOW_STRING(refinement);
	SHOW(tmax);
	return true;
}
