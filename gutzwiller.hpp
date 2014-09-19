/*
 * gutzwiller.hpp
 *
 *  Created on: Sep 11, 2014
 *      Author: Abuenameh
 */

#ifndef GUTZWILLER_HPP_
#define GUTZWILLER_HPP_

#include "cuda_complex.hpp"
#include "configuration.h"

#define L 5
#define nmax 5
#define dim (nmax+1)

__host__ __device__ inline int mod(int i) {
	return (i + L) % L;
}

template<class T>
__host__ __device__ inline T g(int n, int m) {
	return sqrt((T)1.0*(n + 1) * m);
}

__host__ __device__ inline double eps(real* U, int i, int j, int n, int m) {
	return n * U[i] - (m - 1) * U[j];
}

template<class T>
struct parameters {
	T* U;
	T* J;
	T mu;
    T theta;
    complex<T>* f;
    real costh;
    real sinth;
    real cos2th;
    real sin2th;
//	real* U;
//	real mu;
//	real* J;
//    real theta;
//    real costh;
//    real sinth;
//    real cos2th;
//    real sin2th;
};

template<class T>
class Energy {
public:
	__host__ __device__ T operator()(const T *f, unsigned int n,
			void *f_data) const;
};

//double f_nelderMead(unsigned int n, const double *x, double *grad,
//		void *f_data);


#endif /* GUTZWILLER_HPP_ */
