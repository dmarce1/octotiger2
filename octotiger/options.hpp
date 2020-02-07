#pragma once
#include <string>

class options {
public:
	bool global_time;
	double cfl;
	double fgamma;
	double grid_size;
	std::string config_file;
	std::string problem;

	template<class Arc>
	void serialize(Arc &arc, unsigned) {
		arc & global_time;
		arc & cfl;
		arc & fgamma;
		arc & grid_size;
		arc & config_file;
		arc & problem;
	}
	static options global;
	static options& get();
	static void set(options);
	bool process_options(int argc, char *argv[]);
};
