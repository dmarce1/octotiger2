/*
 * sub_array.hpp
 *
 *  Created on: Feb 21, 2020
 *      Author: dmarce1
 */

#ifndef OCTOTIGER_SUB_ARRAY_HPP_
#define OCTOTIGER_SUB_ARRAY_HPP_

#include <octotiger/volume.hpp>

#include <cassert>

using index_type = general_vect<int,NDIM>;


template<class T>
class sub_array {
	std::vector<T> data_;
	volume<int> volume_;
	void initialize() {
		data_.resize(volume_.size());
	}
public:
	sub_array(const volume<int> volume) :
			volume_(volume) {
		initialize();
	}

	const T& operator[](const index_type &I) const {
		assert(volume_.contains(I));
		const auto i = volume_.index(I);
		return data_[i];
	}

	T& operator[](const index_type &I) {
		assert(volume_.contains(I));
		const auto i = volume_.index(I);
		return data_[i];
	}
};


#endif /* OCTOTIGER_SUB_ARRAY_HPP_ */
