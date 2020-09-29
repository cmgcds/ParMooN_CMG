#include "convolution.h"

// Computes discrete version of w(x) = \int u(x-y) * v(y) dy on h-refined grid (a, b]
// Here is assumed that a and b are pieacewise linear, and c is also projected onto space 
// of piecewise linear functions
void linear_convolution(double* u, double* v, double* w, int start_n, int L, double a, double b) {
	int n = (L + 1) * start_n / 2.0 + 1, i, j;
	double h = (b - a) / start_n;
	double* x = allocate_vector(n);
	fill_hrefined_grid(a, b, L, start_n, x);

	for(i = 0; i < n - 1; i++)
		for(j = 0; j <= i; j++)
			w[i] = u[i - j] * v[j] * (x[j + 1] - x[j]);

/*	double* coeff_u0 = allocate_vector(n - 1);
	double* coeff_u1 = allocate_vector(n - 1);
	double* coeff_v0 = allocate_vector(n - 1);
	double* coeff_v1 = allocate_vector(n - 1);
	for(int i = 0; i < n - 1; i++) {
		get_Legendre_coefficients(u + i, u + i + 1, x + i, x + i + 1, coeff_u0 + i, coeff_u1 + i);
		get_Legendre_coefficients(v + i, v + i + 1, x + i, x + i + 1, coeff_v0 + i, coeff_v1 + i);
	} */

}

// Gets Legendre coefficients for linear functions
void get_Legendre_coefficients(double* fi, double* fiplus1, double* xi, double* xiplus1, double* coeff0, double* coeff1) {

}