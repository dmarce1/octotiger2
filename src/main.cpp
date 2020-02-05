
#include <hpx/hpx_init.hpp>

int hpx_main(int argc, char *argv[]) {

	return hpx::finalize();
}

int main(int argc, char *argv[]) {
	std::vector<std::string> cfg = { "hpx.commandline.allow_unknown=1" };
	hpx::init(argc, argv, cfg);
}
