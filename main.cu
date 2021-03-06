/*
 * main.cpp
 *
 *  Created on: Sep 11, 2014
 *      Author: Abuenameh
 */

#include <ctime>

#include "cusimann.cuh"
#include "nelderMead.h"
#include "gutzwiller.hpp"

double f_nelderMead(unsigned int n, const double *x, double *grad,
		void *f_data) {
	return Energy<double>()(x, n, f_data);
}

int main(int argc, char** argv) {
	time_t start = time(NULL);

	real T_0 = 1000, T_min = 0.1;
	const unsigned int n = 2 * L * dim, N = 10;
	const real rho = 0.99;
	size_t sizeFD = n * sizeof(real);
	real *lb, *ub, *cusimann_minimum = (real*) malloc(sizeFD),
			f_cusimann_minimum;
	lb = (real*) malloc(sizeFD);
	unsigned int i;
	for (i = 0; i < n; i++)
		lb[i] = -1;
	ub = (real*) malloc(sizeFD);
	for (i = 0; i < n; i++)
		ub[i] = 1;

	unsigned int n_threads_per_block = 128;//512;//256;
	unsigned int n_blocks = 64;

	double U[L], J[L];
	for (int i = 0; i < L; i++) {
		U[i] = 1;
		J[i] = 0.001;
	}

	real h_d_U[L], h_d_J[L];
	for (int i = 0; i < L; i++) {
		h_d_U[i] = (real)U[i];
		h_d_J[i] = (real)J[i];
	}

	parameters<real> h_d_parms;
	parameters<real>* d_parms;
	real* d_U;
	real* d_J;
	checkCudaErrors(cudaMalloc(&d_U, L*sizeof(real)));
	checkCudaErrors(cudaMemcpy(d_U, h_d_U, L*sizeof(real), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMalloc(&d_J, L*sizeof(real)));
	checkCudaErrors(cudaMemcpy(d_J, h_d_J, L*sizeof(real), cudaMemcpyHostToDevice));

	real theta = 0;

	h_d_parms.U = d_U;
	h_d_parms.J = d_J;
	h_d_parms.mu = 0.5;
	h_d_parms.theta = theta;
	h_d_parms.costh = cos(theta);
	h_d_parms.sinth = sin(theta);
	h_d_parms.cos2th = cos(2*theta);
	h_d_parms.sin2th = sin(2*theta);
	checkCudaErrors(cudaMalloc(&d_parms, sizeof(parameters<real>)));
	checkCudaErrors(
			cudaMemcpy(d_parms, &h_d_parms, sizeof(parameters<real>),
					cudaMemcpyHostToDevice));

	cusimann_optimize(n_threads_per_block, n_blocks, T_0, T_min, N, rho, n, lb,
			ub, Energy<real>(), d_parms, cusimann_minimum, &f_cusimann_minimum);

	printf("cusimann_minimum = [");
	for (i = 0; i < n; i++)
		printf(" %f", cusimann_minimum[i]);
	printf(" ]\n");
	printf("f(cusimann_minimum) = %lf\n", f_cusimann_minimum);

	parameters<double> parms;
	parms.U = U;
	parms.J = J;
	parms.mu = 0.5;
	parms.theta = theta;

	double f_nelderMead_minimum;
	double *nelderMead_minimum = (double*) malloc(n * sizeof(double));
	nelderMead_optimize(n, lb, ub, cusimann_minimum, f_nelderMead, &parms,
			nelderMead_minimum, &f_nelderMead_minimum);

	printf("nelderMead_minimum = [");
	for (i = 0; i < n; i++)
		printf(" %f", nelderMead_minimum[i]);
	printf(" ]\n");
	printf("f(nelderMead_minimum) = %lf\n", f_nelderMead_minimum);

	free(lb);
	free(ub);
	free(cusimann_minimum);
	free(nelderMead_minimum);

	time_t end = time(NULL);

	printf("Runtime: %ld s\n", end-start);

	return 0;
}

