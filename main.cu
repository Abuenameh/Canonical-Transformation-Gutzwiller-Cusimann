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

	real T_0 = 1000, T_min = 0.01;
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

	real U[L], J[L];
	for (int i = 0; i < L; i++) {
		U[i] = 1;
		J[i] = 0.001;
	}
	parameters parms;
	parms.U = U;
	parms.J = J;
	parms.mu = 0.5;

	parameters* d_parms;
	real* d_U;
	real* d_J;
	checkCudaErrors(cudaMalloc(&d_U, L*sizeof(real)));
	checkCudaErrors(cudaMemcpy(d_U, U, L*sizeof(real), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMalloc(&d_J, L*sizeof(real)));
	checkCudaErrors(cudaMemcpy(d_J, J, L*sizeof(real), cudaMemcpyHostToDevice));

	real theta = 0;

	parms.U = d_U;
	parms.J = d_J;
	parms.mu = 0.5;
	parms.theta = theta;
	parms.costh = cos(theta);
	parms.sinth = sin(theta);
	parms.cos2th = cos(2*theta);
	parms.sin2th = sin(2*theta);
	checkCudaErrors(cudaMalloc(&d_parms, sizeof(parameters)));
	checkCudaErrors(
			cudaMemcpy(d_parms, &parms, sizeof(parameters),
					cudaMemcpyHostToDevice));

	printf("Optimize\n");
	cusimann_optimize(n_threads_per_block, n_blocks, T_0, T_min, N, rho, n, lb,
			ub, Energy<real>(), d_parms, cusimann_minimum, &f_cusimann_minimum);

	printf("cusimann_minimum = [");
	for (i = 0; i < n; i++)
		printf(" %f", cusimann_minimum[i]);
	printf(" ]\n");
	printf("f(cusimann_minimum) = %lf\n", f_cusimann_minimum);

	parms.U = U;
	parms.J = J;

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

template<class T>
class Ackley {
public:
	__host__ __device__ T operator()(const T *x, unsigned int n,
			void *f_data) const {
		T f_x = 0.0f;
		T aux = exp(-0.2f);

		int i;
		for (i = 0; i < n - 1; i++)
			f_x = f_x + aux * sqrt(pow(x[i], 2) + pow(x[i + 1], 2))
					+ 3.0f * (cos(2.0f * x[i]) + sin(2.0f * x[i + 1]));

		printf("Iter\n");
		return f_x;
	}
};

double f_nelderMead2(unsigned int n, const double *x, double *grad,
		void *f_data) {
	return Ackley<double>()(x, n, f_data);
}

int main2() {
	real T_0 = 1000, T_min = 0.1;
	const unsigned int n = 5, N = 100;
	const real rho = 0.99;
	size_t sizeFD = n * sizeof(real);
	real *lb, *ub, *cusimann_minimum = (real*) malloc(sizeFD),
			f_cusimann_minimum;
	lb = (real*) malloc(sizeFD);
	unsigned int i;
	for (i = 0; i < n; i++)
		lb[i] = -30;
	ub = (real*) malloc(sizeFD);
	for (i = 0; i < n; i++)
		ub[i] = 30;

	unsigned int n_threads_per_block = 256;
	unsigned int n_blocks = 64;

	cusimann_optimize(n_threads_per_block, n_blocks, T_0, T_min, N, rho, n, lb,
			ub, Ackley<real>(), NULL, cusimann_minimum, &f_cusimann_minimum);

	cudaDeviceSynchronize();
	printf("cusimann_minimum = [");
	for (i = 0; i < n; i++)
		printf(" %f", cusimann_minimum[i]);
	printf(" ]\n");
	printf("f(cusimann_minimum) = %lf\n", f_cusimann_minimum);

	double f_nelderMead_minimum;
	double *nelderMead_minimum = (double*) malloc(n * sizeof(double));
//	nelderMead_optimize(n, lb, ub, cusimann_minimum, f_nelderMead2, NULL, nelderMead_minimum, &f_nelderMead_minimum);

	printf("nelderMead_minimum = [");
	for (i = 0; i < n; i++)
		printf(" %f", nelderMead_minimum[i]);
	printf(" ]\n");


	free(lb);
	free(ub);
	free(cusimann_minimum);
	free(nelderMead_minimum);

	return EXIT_SUCCESS;
}
