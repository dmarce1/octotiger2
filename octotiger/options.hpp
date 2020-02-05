#pragma once
#include <string>

class options {
public:
	std::string config_file;
	bool global_time;
	double fgamma;
	double grid_size;
	double cfl;

	template<class Arc>
	void serialize(Arc &arc, unsigned) {
		arc & cfl;
		arc & config_file;
		arc & global_time;
		arc & fgamma;
		arc & grid_size;
	}
	static options global;
	static options& get();
	static void set(options);
	bool process_options(int argc, char *argv[]);
};
