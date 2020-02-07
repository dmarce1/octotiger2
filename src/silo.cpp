/*
 * silo.cpp
 *
 *  Created on: Feb 6, 2020
 *      Author: dmarce1
 */
#include <octotiger/silo.hpp>
#include <octotiger/fixed_real.hpp>

#include <hpx/lcos/local/mutex.hpp>

#include <unordered_map>

using node_type = general_vect<fixed_real, NDIM>;
using zone_map_type = std::pair<full_state,std::array<int,NCHILD>>;

struct node_hash {
	std::size_t operator()(const node_type &n) const {
		std::size_t key = 0;
		for (int dim = 0; dim < NDIM; dim++) {
			key ^= n[dim].get_int();
		}
		return key;
	}
};

static hpx::lcos::local::mutex mtx_;
static std::unordered_map<node_type, int, node_hash> node_map_;
static int next_node_num_;
static std::vector<zone_map_type> zone_map_;

void silo_begin() {
	zone_map_.clear();
	node_map_.clear();
	next_node_num_ = 0;
}

void silo_add_zones(const std::vector<silo_zone> &zones) {
	std::lock_guard<hpx::lcos::local::mutex> lock(mtx_);

	for (const auto &z : zones) {
		std::array<int, NCHILD> indices;
		for (int ci = 0; ci < NCHILD; ci++) {
			const auto iter = node_map_.find(z.nodes[ci]);
			int i;
			if (iter == node_map_.end()) {
				i = next_node_num_++;
				node_map_[z.nodes[ci]] = i;
			} else {
				i = iter->second;
			}
			indices[ci] = i;
		}
		zone_map_type zone;
		zone.first = z.state;
		zone.second = indices;
		zone_map_.push_back(zone);
	}

}

void silo_end(const std::string &fname) {
	DBfile *db = DBCreateReal(fname.c_str(), DB_CLOBBER, DB_LOCAL, "Octo-Tiger II", DB_HDF5);

	const int nnodes = node_map_.size();
	const int nzones = zone_map_.size();
	std::vector<int> zones;
	double *coords[NDIM];
	char *coordnames[NDIM];

	for (int dim = 0; dim < NDIM; dim++) {
		coords[dim] = new double[nnodes];
		coordnames[dim] = new char[2];
		coordnames[dim][0] = 'x' + dim;
		coordnames[dim][1] = '\0';
	}

	for (auto iter = node_map_.begin(); iter != node_map_.end(); iter++) {
		for (int dim = 0; dim < NDIM; dim++) {
			coords[dim][iter->second] = iter->first[dim];
		}
	}
#if( NDIM == 1)
		const int shapetypes[1] = {DB_ZONETYPE_BEAM};
#elif( NDIM ==2)
	const int shapetypes[1] = { DB_ZONETYPE_QUAD };
#elif(NDIM==1)
		const int shapetypes[1] = { DB_ZONETYPE_HEX };
#else
#error
#endif
	const int shapesize[1] = { NDIM + 1 };
	const int shapecount[1] = { nzones };

	for (const auto &z : zone_map_) {
#if(NDIM==3)
		zones.push_back(z.second[0]);
		zones.push_back(z.second[1]);
		zones.push_back(z.second[3]);
		zones.push_back(z.second[2]);
		zones.push_back(z.second[4]);
		zones.push_back(z.second[5]);
		zones.push_back(z.second[7]);
		zones.push_back(z.second[6]);
#elif(NDIM==2)
		zones.push_back(z.second[0]);
		zones.push_back(z.second[1]);
		zones.push_back(z.second[3]);
		zones.push_back(z.second[2]);
#elif(NDIM==1)
		zones.push_back(z.second[0]);
		zones.push_back(z.second[1]);
#else
#error
#endif
	}
	DBPutZonelist2(db, "zonelist", nzones, NDIM, zones.data(), zones.size(), 0, 0, 0, shapetypes, shapesize, shapecount, 1, NULL);
	DBPutUcdmesh(db, "mesh", NDIM, coordnames, coords, nnodes, nzones, "zonelist", NULL, DB_DOUBLE, NULL);





	for (int dim = 0; dim < NDIM; dim++) {
		delete[] coordnames[dim];
		delete[] coords[dim];
	}
	for (int dim = 0; dim < NDIM; dim++) {
		delete[] coords[dim];
	}
	DBClose(db);
}
