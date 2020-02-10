#pragma once
#include <string>

class options {
public:
	bool global_time;
	double cfl;
	double fgamma;
	double output_freq;
	double tmax;
	int grid_size;
	int max_level;
	std::string config_file;
	std::string problem;
	std::string refinement;

	template<class Arc>
	void serialize(Arc &arc, unsigned) {
		arc & global_time;
		arc & cfl;
		arc & fgamma;
		arc & output_freq;
		arc & tmax;
		arc & grid_size;
		arc & max_level;
		arc & config_file;
		arc & problem;
		arc & refinement;
	}
	static options global;
	static options& get();
	static void set(options);
	bool process_options(int argc, char *argv[]);
};
