#ifndef MPIMIS_CONVOLUTION
#define MPIMIS_CONVOLUTION

#include "basictools.h"

// Computes discrete version of w(x) = \int u(x-y) * v(y) dy on h-refined grid (a, b]
// Here is assumed that a and b are pieacewise linear, and c is also projected onto space 
// of piecewise linear functions
void linear_convolution(double* u, double* v, double* w, int start_n, int L, double a, double b);

// Gets Legendre coefficients for linear functions
void get_Legendre_coefficients(double* fi, double* fiplus1, double* xi, double* xiplus1, double* coeff0, double* coeff1);

#endif // MPIMIS_CONVOLUTION

