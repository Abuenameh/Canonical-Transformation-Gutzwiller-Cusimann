/*
 * energy.cu
 *
 *  Created on: Sep 18, 2014
 *      Author: Abuenameh
 */

#include "gutzwiller.hpp"

template<class T>
__host__ __device__ T Energy<T>::operator ()(const T *x, unsigned int n,
	void *f_data) const {
	parameters* parms = (parameters*) f_data;
	real* U = parms->U;
	real* J = parms->J;
	real mu = parms->mu;
	real costh = parms->costh;
	real sinth = parms->sinth;

	typedef typename complextype<T>::type complex_t;

	complex_t expth = complextype<T>::make_complex(costh, sinth);
	complex_t expmth = ~expth;
	complex_t exp2th = expth * expth;
	complex_t expm2th = ~exp2th;

	complex_t Ec = complex_t::zero();

	const complex_t * f[L];
	T norm2[L];
	for (int i = 0; i < L; i++) {
		f[i] = reinterpret_cast<const complex_t*>(&x[2 * i * dim]);
		norm2[i] = 0;
		for (int n = 0; n <= nmax; n++) {
			norm2[i] += norm(f[i][n]);
		}
	}

	for (int i = 0; i < L; i++) {

		int k1 = mod(i - 2);
		int j1 = mod(i - 1);
		int j2 = mod(i + 1);
		int k2 = mod(i + 2);

		complex_t E0 = complex_t::zero();
		complex_t E1j1 = complex_t::zero();
		complex_t E1j2 = complex_t::zero();
		complex_t E2j1 = complex_t::zero();
		complex_t E2j2 = complex_t::zero();
		complex_t E3j1 = complex_t::zero();
		complex_t E3j2 = complex_t::zero();
		complex_t E4j1j2 = complex_t::zero();
		complex_t E4j1k1 = complex_t::zero();
		complex_t E4j2k2 = complex_t::zero();
		complex_t E5j1j2 = complex_t::zero();
		complex_t E5j1k1 = complex_t::zero();
		complex_t E5j2k2 = complex_t::zero();

        for (int n = 0; n <= nmax; n++) {
            E0 += (0.5 * U[i] * n * (n - 1) - mu * n) * ~f[i][n] * f[i][n];

            if (n < nmax) {
                E1j1 += -J[j1] * expth * g(n, n + 1) * ~f[i][n + 1] * ~f[j1][n]
                        * f[i][n] * f[j1][n + 1];
                E1j2 += -J[i] * expmth * g(n, n + 1) * ~f[i][n + 1] * ~f[j2][n] * f[i][n]
                        * f[j2][n + 1];

            }

        }

		Ec += E0 / norm2[i];

		Ec += E1j1 / (norm2[i] * norm2[j1]);
		Ec += E1j2 / (norm2[i] * norm2[j2]);

		Ec += E2j1 / (norm2[i] * norm2[j1]);
		Ec += E2j2 / (norm2[i] * norm2[j2]);

		Ec += E3j1 / (norm2[i] * norm2[j1]);
		Ec += E3j2 / (norm2[i] * norm2[j2]);

		Ec += E4j1j2 / (norm2[i] * norm2[j1] * norm2[j2]);
		Ec += E4j1k1 / (norm2[i] * norm2[j1] * norm2[k1]);
		Ec += E4j2k2 / (norm2[i] * norm2[j2] * norm2[k2]);

		Ec += E5j1j2 / (norm2[i] * norm2[j1] * norm2[j2]);
		Ec += E5j1k1 / (norm2[i] * norm2[j1] * norm2[k1]);
		Ec += E5j2k2 / (norm2[i] * norm2[j2] * norm2[k2]);
	}

	return Ec.real();
}

template class Energy<float> ;
template class Energy<double> ;

