/*
 * complex.hpp
 *
 *  Created on: Sep 11, 2014
 *      Author: Abuenameh
 */

#ifndef COMPLEX_HPP_
#define COMPLEX_HPP_

//#include <complex>
//using std::complex;

#include "cudacomplex.cuh"

//#ifdef __CUDACC__
//#else
//#endif

template<class T>
struct complextype {
	typedef T type;
	__host__ __device__ static inline type make_complex(T a, T b);
};

template<>
struct complextype<float> {
	typedef singlecomplex type;
	__host__ __device__ static inline type make_complex(float a, float b) {
		return make_singlecomplex(a, b);
	}
};

template<>
struct complextype<double> {
	typedef doublecomplex type;
	__host__ __device__ static inline type make_complex(double a, double b) {
		return make_doublecomplex(a, b);
	}
};

//template<>
//struct complextype<float> {
//#ifdef __CUDACC__
//	typedef singlecomplex type;
//	static inline type make_complex(float a, float b) {
//		return make_singlecomplex(a, b);
//	}
//#else
//	typedef complex<float> type;
//	static inline type make_complex(float a, float b) {
//		return type(a, b);
//	}
//#endif
//};

//template<>
//struct complextype<double> {
//#ifdef __CUDACC__
//	typedef doublecomplex type;
//	static inline type make_complex(double a, double b) {
//		return make_doublecomplex(a, b);
//	}
//#else
//	typedef complex<double> type;
//	static inline type make_complex(double a, double b) {
//		return type(a, b);
//	}
//#endif
//};

//template<class T>
//complex<T> operator~(const complex<T> a) {
//	return conj(a);
//}




#endif /* COMPLEX_HPP_ */
