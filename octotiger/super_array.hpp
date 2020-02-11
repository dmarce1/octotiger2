/*
 * super_array.hpp
 *
 *  Created on: Feb 5, 2020
 *      Author: dmarce1
 */

#ifndef OCTOTIGER_SUPER_ARRAY_HPP_
#define OCTOTIGER_SUPER_ARRAY_HPP_

#include <octotiger/volume.hpp>

#include <hpx/lcos/local/shared_mutex.hpp>

#include <cassert>
#include <set>

using index_type = general_vect<int,NDIM>;

template<class T>
class super_array {
	std::vector<T> data_;
	volume<int> volume_;
	std::set<volume<int>> sub_volumes_;
	mutable hpx::lcos::local::shared_mutex mtx_;

	void resize(const volume<int> &new_volume) {
		if (new_volume != volume_) {
			std::vector<T> new_data(new_volume.size());
			for (auto I0 = volume_.begin(); I0 != volume_.end(); volume_.inc_index(I0)) {
				auto I1 = I0 - volume_.begin() + new_volume.begin();
				new_data[new_volume.index(I1)] = data_[volume_.index(I0)];
			}
			data_ = std::move(new_data);
			volume_ = std::move(new_volume);
		}
	}

public:

	hpx::lcos::local::shared_mutex& get_mutex() const {
		return mtx_;
	}

	void add_volume(const volume<int> &this_vol) {
		std::lock_guard<hpx::lcos::local::shared_mutex> lock(mtx_);
		if (!data_.empty()) {
			const auto new_volume = volume_.union_(this_vol);
			resize(new_volume);
		} else {
			volume_ = this_vol;
			data_.resize(volume_.size());
		}
		sub_volumes_.insert(this_vol);
	}

	void remove_volume(const volume<int> &this_vol) {
		std::lock_guard<hpx::lcos::local::shared_mutex> lock(mtx_);
		assert(sub_volumes_.find(this_vol) != sub_volumes_.end());
		sub_volumes_.erase(this_vol);
	}

	void shrink_to_fit() {
		std::lock_guard<hpx::lcos::local::shared_mutex> lock(mtx_);
		volume<int> new_volume;
		for (const auto &sv : sub_volumes_) {
			new_volume = new_volume.union_(sv);
		}
		resize(new_volume);
	}

	const T& operator[](const index_type& I) const {
		return data_[volume_.index(I)];
	}

	T& operator[](const index_type& I) {
		return data_[volume_.index(I)];
	}

};

#endif /* OCTOTIGER_SUPER_ARRAY_HPP_ */

