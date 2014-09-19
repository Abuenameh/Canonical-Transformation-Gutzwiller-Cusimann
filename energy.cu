/*
 * energy.cu
 *
 *  Created on: Sep 18, 2014
 *      Author: Abuenameh
 */

#include <stdio.h>

#include "gutzwiller.hpp"
#include "cuda_complex.hpp"

template<class T>
__host__ __device__ T Energy<T>::operator ()(const T *x, unsigned int n,
	void *f_data) const {

	parameters<T>* parms = (parameters<T>*) f_data;
	T* U = parms->U;
	T* J = parms->J;
	T mu = parms->mu;
	T theta = parms->theta;
//	T costh = parms->costh;
//	T sinth = parms->sinth;

//	parameters* parms = (parameters*) f_data;
//	real* U = parms->U;
//	real* J = parms->J;
//	real mu = parms->mu;
//	real theta = parms->theta;
//	real costh = parms->costh;
//	real sinth = parms->sinth;

	complex<T> expth = exp(complex<T>(0, 1) * theta);
	complex<T> expmth = ~expth;
	complex<T> exp2th = expth * expth;
	complex<T> expm2th = ~exp2th;

	complex<T> E = 0;
	T Er = 0;
	T Eri[L];

	const complex<T> * f[L];
	T norm2[L];
	complex<T> f2[dim];
	for (int n = 0; n <= nmax; n++) {
		f2[n] = complex<T>(x[2*n],x[2*n+1]);
	}
	for (int i = 0; i < L; i++) {
		f[i] = reinterpret_cast<const complex<T>*>(&x[2 * i * dim]);
//		f[i] = f2;
		norm2[i] = 0;
		for (int n = 0; n <= nmax; n++) {
			norm2[i] += norm(f[i][n]);
		}
	}

//	typedef typename complextype<T>::type complex_t;
//
//	complex_t expth = complextype<T>::make_complex(costh, sinth);
//	complex_t expmth = ~expth;
//	complex_t exp2th = expth * expth;
//	complex_t expm2th = ~exp2th;
//
//	complex_t Ec = complex_t::zero();

//	const complex_t * f[L];
//	T norm2[L];
//	for (int i = 0; i < L; i++) {
//		f[i] = reinterpret_cast<const complex_t*>(&x[2 * i * dim]);
//		norm2[i] = 0;
//		for (int n = 0; n <= nmax; n++) {
//			norm2[i] += norm(f[i][n]);
//		}
//	}

	for (int i = 0; i < L; i++) {

		int k1 = mod(i - 2);
		int j1 = mod(i - 1);
		int j2 = mod(i + 1);
		int k2 = mod(i + 2);

		complex<T> E0 = 0;
		complex<T> E1j1 = 0;
		complex<T> E1j2 = 0;
		complex<T> E2j1 = 0;
		complex<T> E2j2 = 0;
		complex<T> E3j1 = 0;
		complex<T> E3j2 = 0;
		complex<T> E4j1j2 = 0;
		complex<T> E4j1k1 = 0;
		complex<T> E4j2k2 = 0;
		complex<T> E5j1j2 = 0;
		complex<T> E5j1k1 = 0;
		complex<T> E5j2k2 = 0;

//		complex_t E0 = complex_t::zero();
//		complex_t E1j1 = complex_t::zero();
//		complex_t E1j2 = complex_t::zero();
//		complex_t E2j1 = complex_t::zero();
//		complex_t E2j2 = complex_t::zero();
//		complex_t E3j1 = complex_t::zero();
//		complex_t E3j2 = complex_t::zero();
//		complex_t E4j1j2 = complex_t::zero();
//		complex_t E4j1k1 = complex_t::zero();
//		complex_t E4j2k2 = complex_t::zero();
//		complex_t E5j1j2 = complex_t::zero();
//		complex_t E5j1k1 = complex_t::zero();
//		complex_t E5j2k2 = complex_t::zero();

		complex<T> Eg;
        for (int n = 0; n <= nmax; n++) {
////        	E0 += f[i][n];
//        	            complex<T> Ee = /*((T)0.5 * U[i] * n * (n - 1) - mu * n) */ ~f[i][n] * ~f[i][n];
//        	            complex<T> Ef = Ee;///*((T)0.5 * U[i] * n * (n - 1) - mu * n) */ ~f[i][n] * ~f[i][n];
//        	            Eg += Ee;
//        	printf("%d %d\n", i, n);
//        	complex<T> asd = f[i][n];
//        	const complex<T>* poi = f[0];
//        	const complex<T>* lkj = f[i];
//        	complex<T> oiu = poi[n];
//        	complex<T> kjh = lkj[n];
//        	complex<T> mnb = f[i][n];
//        	complex<T> zxc = f[L-1][nmax];
//        	complex<T> sdf = ~asd*asd;
//        	complex<T> wer = ~f[i][n]*f[i][n];
//            E0 = ((T)0.5 * U[i] * n * (n - 1) - mu * n) * ~f[i][n] * f[i][n];
//            E0 += ((T)0.5 * U[i] * n * (n - 1) - mu * n) * ~f2[n] * f2[n];
            E0 = ((T)0.5 * U[i] * n * (n - 1) - mu * n) * ~f[i][n] * f[i][n];

            if (n < nmax) {
//            	complex_t qwe = -J[j1]*expth*g(n,n+1)*f[i][n+1];//*~f[j1][n]*f[i][n]*f[j1][n+1];
//            	E1j1 += qwe;
//            	E1j1 += -J[j1]*expth*g(n,n+1)*f[i][n+1];
//            	E1j1 += -J[j1]*expth*g(n,n+1);
//            	E1j1 += expth*f[i][n+1];
//                E1j1 += -J[j1] * expth * g<T>(n, n + 1) * ~f[i][n + 1] * ~f[j1][n]
//                        * f[i][n] * f[j1][n + 1];
//                E1j2 += -J[i] * expmth * g<T>(n, n + 1) * ~f[i][n + 1] * ~f[j2][n] * f[i][n]
//                        * f[j2][n + 1];
            }

        }

		E += E0 / norm2[i];
//		Eri[i] = E0.real() / norm2[i] + E1j1.real() / (norm2[i] * norm2[j1]) + E1j2.real() / (norm2[i] * norm2[j2]);

		E += E1j1 / (norm2[i] * norm2[j1]);
		E += E1j2 / (norm2[i] * norm2[j2]);
//				Er += E1j1.real() / (norm2[i] * norm2[j1]);
//				Er += E1j2.real() / (norm2[i] * norm2[j2]);

//		E += E2j1 / (norm2[i] * norm2[j1]);
//		E += E2j2 / (norm2[i] * norm2[j2]);
//
//		E += E3j1 / (norm2[i] * norm2[j1]);
//		E += E3j2 / (norm2[i] * norm2[j2]);
//
//		E += E4j1j2 / (norm2[i] * norm2[j1] * norm2[j2]);
//		E += E4j1k1 / (norm2[i] * norm2[j1] * norm2[k1]);
//		E += E4j2k2 / (norm2[i] * norm2[j2] * norm2[k2]);
//
//		E += E5j1j2 / (norm2[i] * norm2[j1] * norm2[j2]);
//		E += E5j1k1 / (norm2[i] * norm2[j1] * norm2[k1]);
//		E += E5j2k2 / (norm2[i] * norm2[j2] * norm2[k2]);
	}
//    printf("Here\n");

	return E.real();
//	return Eri[0];
//	return 0;
}

template class Energy<float> ;
template class Energy<double> ;

