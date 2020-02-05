/*
 * volume.hpp
 *
 *  Created on: Feb 1, 2020
 *      Author: dmarce1
 */

#ifndef OCTOTIGER_VOLUME_HPP_
#define OCTOTIGER_VOLUME_HPP_

#include <octotiger/vect.hpp>

#include <cassert>

template<class T>
class volume {
	general_vect<T, NDIM> begin_;
	general_vect<T, NDIM> end_;

public:
	template<class Arc>
	inline void serialize(Arc &a, unsigned) {
		a & begin_;
		a & end_;
	}

	volume() {
		for (int dim = 0; dim < NDIM; dim++) {
			begin(dim) = end(dim) = 0;
		}
	}

	bool empty() const {
		bool rc = false;
		for (int dim = 0; dim < NDIM; dim++) {
			if (end(dim) <= begin(dim)) {
				rc = true;
				break;
			}
		}
		return rc;
	}

	volume(const std::array<T, NDIM> &b, const std::array<T, NDIM> &e) {
		for (int dim = 0; dim < NDIM; dim++) {
			begin(dim) = b[dim];
			end(dim) = e[dim];
		}
	}

	inline T begin(int dim) const {
		return begin_[dim];
	}

	inline T& begin(int dim) {
		return begin_[dim];
	}

	inline T end(int dim) const {
		return end_[dim];
	}

	inline T& end(int dim) {
		return end_[dim];
	}

	volume<T> union_(const volume<T> &other) const {
		volume<T> u;
		if (empty()) {
			u = other;
		} else if (other.empty()) {
			u = *this;
		} else {
			using std::min;
			using std::max;
			for (int dim = 0; dim < NDIM; dim++) {
				u.begin(dim) = min(begin(dim), other.begin(dim));
				u.end(dim) = max(end(dim), other.end(dim));
			}
		}
		return u;
	}

	bool operator<(const volume<T> &other) const {
		bool rc = false;
		for (int dim = 0; dim < NDIM; dim++) {
			if (begin(dim) < other.begin(dim)) {
				rc = true;
				break;
			} else if (begin(dim) > other.begin(dim)) {
				rc = false;
				break;
			}
			if (end(dim) < other.end(dim)) {
				rc = true;
				break;
			} else if (end(dim) > other.end(dim)) {
				rc = false;
				break;
			}
		}
		return rc;
	}

	bool operator==(const volume<T> &other) const {
		bool rc = true;
		for (int dim = 0; dim < NDIM; dim++) {
			if (begin(dim) != other.begin(dim)) {
				rc = false;
				break;
			}
			if (end(dim) != other.end(dim)) {
				rc = false;
				break;
			}
		}
		return rc;
	}

	inline bool operator!=(const volume<T> &other) const {
		return !(*this == other);
	}

	general_vect<T, NDIM> begin() const {
		return begin_;
	}

	general_vect<T, NDIM> end() const {
		return end_;
	}

	volume<T> expand(const T& d) {
		volume<T> v;
		for( int dim = 0; dim < NDIM; dim++) {
			v.begin(dim) = begin(dim) - d;
			v.end(dim) = end(dim) + d;
		}
		return v;
	}

	T index(const general_vect<T, NDIM> &I) const {
		T i = I[0] - begin(0);
		for (int dim = 1; dim < NDIM; dim++) {
			i *= end(dim) - begin(dim);
			i += I[dim] - begin(dim);
		}
		return i;
	}

	void inc_index(general_vect<T, NDIM> &I) {
		int dim = 0;
		while (++I[dim] == end(dim)) {
			if (dim != NDIM - 1) {
				I[dim] = 0;
			} else {
				I = end();
			}
			dim++;
		}
	}

	T size() const {
		T s = end(0) - begin(0);
		for (int dim = 1; dim < NDIM; dim++) {
			s *= end(dim) - begin(dim);
		}
		return s;
	}

	std::array<T, NDIM> sizes() const {
		std::array<T, NDIM> s;
		for (int dim = 0; dim < NDIM; dim++) {
			s[dim] = end(dim) - begin(dim);
		}
		return s;
	}

	std::array<T, NDIM> strides() const {
		std::array<T, NDIM> s;
		s[NDIM - 1] = 1;
		for (int dim = NDIM - 2; dim >= 0; dim--) {
			s[dim] = s[dim + 1] * (begin(dim + 1) - end(dim + 1));
		}
		return s;
	}

};

#endif /* OCTOTIGER_VOLUME_HPP_ */
