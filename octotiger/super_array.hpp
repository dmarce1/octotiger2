/*
 * super_array.hpp
 *
 *  Created on: Feb 5, 2020
 *      Author: dmarce1
 */

#ifndef OCTOTIGER_SUPER_ARRAY_HPP_
#define OCTOTIGER_SUPER_ARRAY_HPP_

#include <octotiger/sub_array.hpp>

#include <hpx/lcos/local/mutex.hpp>

#include <set>

template<class T>
class super_array {
	std::vector<T> data_;
	volume<int> volume_;
	std::set<volume<int>> sub_volumes_;
	mutable hpx::lcos::local::mutex mtx_;

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

	void add_volume(const volume<int> &this_vol) {
		std::lock_guard<hpx::lcos::local::mutex> lock(mtx_);
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
		std::lock_guard<hpx::lcos::local::mutex> lock(mtx_);
		assert(sub_volumes_.find(this_vol) != sub_volumes_.end());
		sub_volumes_.erase(this_vol);
	}

	void shrink_to_fit() {
		std::lock_guard<hpx::lcos::local::mutex> lock(mtx_);
		volume<int> new_volume;
		for (const auto &sv : sub_volumes_) {
			new_volume = new_volume.union_(sv);
		}
		resize(new_volume);
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

	void set_subarray(const sub_array<T> &sub) {
		assert(volume_.contains(sub.volume_));
		for (auto I = sub.volume_.begin(); I != sub.volume_.end(); sub.volume_.inc_index(I)) {
			(*this)[I] = sub[I];
		}
	}

	sub_array<T> get_subarray(const volume<int> &vol) {
		sub_array<T> sub(vol);
		for (auto I = vol.begin(); I != vol.end(); vol.inc_index(I)) {
			sub[I] = (*this)[I];
		}
		return sub;
	}

	sub_array<T> get_restricted_subarray(const volume<int> &vol) {
		const auto rvol = vol.half();
		sub_array<T> R(rvol);
		for (auto I = rvol.begin(); I != rvol.end(); rvol.inc_index(I)) {
			R[I] = (*this)[I * 2] / real(NCHILD);
			for (int ci = 1; ci < NCHILD; ci++) {
				auto J = I * 2;
				for (int dim = 0; dim < NDIM; dim++) {
					if (((ci >> dim) & 1) == 1) {
						J[dim]++;
					}
				}
				R[I] = R[I] + (*this)[J] / real(NCHILD);
			}
		}
		return R;
	}

	sub_array<T> get_prolonged_subarray(const volume<int> &pvol) {
		sub_array<T> P(pvol);
		for (auto I = pvol.begin(); I != pvol.end(); pvol.inc_index(I)) {
			const auto J = I / 2;
			P[I] = (*this)[J];
		}
		return P;
	}

};

#endif /* OCTOTIGER_SUPER_ARRAY_HPP_ */

