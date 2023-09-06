

/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"
#include "cblas.h"
#define ALPHA 	1
#define BETA 	1
/*
 * Add your BLAS implementation here
 */
double* my_solver(int N, double *A, double *B) {
	double *AxBxAt = (double *)malloc(N * N * sizeof(double));
	double *BtxBt = (double *)malloc(N * N * sizeof(double));
	if(AxBxAt == NULL || BtxBt == NULL) {
		fprintf(stderr, "Malloc failed!\n");
		return NULL;
	}
	int power_of_N = N * N;
	for(register int i = 0; i < power_of_N; ++i) {
		AxBxAt[i] = B[i];
	}
	cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, N, ALPHA, A, N, AxBxAt, N);
	cblas_dtrmm(CblasRowMajor, CblasRight, CblasUpper, CblasTrans, CblasNonUnit, N, N, ALPHA, A, N, AxBxAt, N);
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, N, N, N, ALPHA, B, N, B, N, BETA, AxBxAt, N);

	free(BtxBt);
	return AxBxAt;
}

