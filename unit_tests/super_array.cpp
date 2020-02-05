#include <hpx/hpx_init.hpp>

#include <octotiger/super_array.hpp>

int hpx_main(int argc, char *argv[]) {

	super_array<real> super;
#if(NDIM==1)
	volume<int> vol1( { { 0 } }, { { 3 } });
	volume<int> vol2( { { 2 } }, { { 5 } });
#elif(NDIM==2)
	volume<int> vol1( { { 0, 0 } }, { { 3, 3 } });
	volume<int> vol2( { { 0, 2 } }, { { 3, 6 } });
#else
	volume<int> vol1( { { 0, 0, 1 } }, { { 3, 3, 4 } });
	volume<int> vol2( { { 0, 3, 1 } }, { { 3, 6, 4 } });
#endif

	super.add_volume(vol1);
	super.add_volume(vol2);

	return hpx::finalize();
}

int main(int argc, char *argv[]) {
	std::vector<std::string> cfg = { "hpx.commandline.allow_unknown=1" };
	hpx::init(argc, argv, cfg);
}
