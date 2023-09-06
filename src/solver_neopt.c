/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"

/*
 * Add your unoptimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
	double *AxBxAt = (double *) calloc(N * N, sizeof(double));
	double *AxB = (double *) calloc(N * N, sizeof(double));
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			for(int k = i; k < N; k++) {
				AxB[i * N + j] += A[i * N + k] * B[k * N + j];
			}
		}
	}
	for(int i = 0; i < N; i ++) {
		for(int j = 0; j < N; j++) {
			for(int k = j; k < N; k ++) {
				AxBxAt[i * N + j] += AxB[i * N + k] * A[j * N + k];
			}
		}
	}
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			for(int k = 0; k < N; k++) {
				AxBxAt[i * N + j] += B[k * N + i] * B[j * N + k];
			}
		}
	}
	free(AxB);
	return AxBxAt;
}
