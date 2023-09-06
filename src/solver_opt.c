/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"

/*
 * Add your optimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
	register double *AxBxAt = (double *) calloc(N * N, sizeof(double));
	register double *AxB = (double *) calloc(N * N, sizeof(double));
	register int i, j, k;
	register int N_cp = N;
	register double * AxB_linie = AxB;
	register double *AxBxAt_linie = AxBxAt;
	register double *AxBxAt_coloana = AxBxAt;
	register double * AxB_ptr2 = AxB;
	register double *B_linie = B;
	register double *B_coloana = B;
	register double *A_linie = A;
	register double *A_ptr2 = A;

	for(i = 0; i < N_cp; ++i) {
		B_coloana = B_linie;
		for(k = i; k < N_cp; ++k) {
			register double value_A = *(A_linie + k);
			AxB_ptr2 = AxB_linie;
	
			for(j = 0; j < N_cp; ++j) {
				(*AxB_ptr2) += (value_A) * (*B_coloana);
				AxB_ptr2 ++;
				B_coloana ++;
			}
		}
		A_linie += N_cp;
		B_linie += N_cp;
		AxB_linie += N_cp;
	}
	AxB_linie = AxB;
	AxB_ptr2 = AxB;
	B_linie = B;
	A_linie = A;
	for(i = 0; i < N_cp; ++i) {
		A_linie = A;
		for(j = 0; j < N_cp; ++j) {
			register double value_AxBxAt = 0;
			AxB_ptr2 = AxB_linie + j;
			A_ptr2 = A_linie + j;

			for(k = j; k < N_cp; ++k ) {
				value_AxBxAt += (*AxB_ptr2) * (*A_ptr2);
				AxB_ptr2 ++;
				A_ptr2 ++;
			}
			(*AxBxAt) = value_AxBxAt;

			AxBxAt += 1;	

			A_linie += N;
		}
		AxB_linie += N_cp;		
	}
	AxBxAt = AxBxAt_linie;

	for(i = 0; i < N_cp; ++i) {
		B_linie = B;
		for(j = 0; j < N_cp; ++j) {
			B_coloana = B + i;
			register double sum = (*AxBxAt_coloana);
			for(k = 0; k < N_cp; ++k) {
				sum += (*B_coloana) * (*B_linie);
				B_linie ++;
				B_coloana += N_cp;
			}
			(*AxBxAt_coloana) = sum;
			AxBxAt_coloana ++;
		}

	}
	free(AxB);
	return AxBxAt;	
}

