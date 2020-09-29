#ifndef RKMATRIX
#define RKMATRIX
#include "basictools.h"
// #include <mkl.h>

#define RK_SVD_EPS 1.0e-16
#define RK_SVD_EPS2 1.0e-32
#define RK_SVD_EPS3 1.0e-64

typedef class rkmatrix* prkmatrix;

/*  An rkmatrix with paramters k,rows,cols is a factorisation
a * transpose(b) with (rows x k) matrix 'a' and (cols x k) matrix 'b'. */

class rkmatrix {
public:
	rkmatrix() { rows = 0; cols = 0; rk = 0; k = 0; a = 0x0; b = 0x0; } ;
	rkmatrix(int k_new, int rows_new, int cols_new) { k = k_new; rk = k_new; rows = rows_new; cols = cols_new; a = allocate_vector(rows * k); b = allocate_vector(cols * k); } ;
	~rkmatrix() { if(a) free_vector(a); if(b) free_vector(b); } ;

	int k; // allocated rank
	int rk; // real rank
	int rows;
	int cols;
	double* a;
	double* b;

	void rkmatrix_times_vector(EvalMode mode, double* b, double* c, double alpha, double beta);
	void rkmatrix_times_matrix(EvalMode mode, double* B, double* C, int n, double alpha, double beta);
	double get_frobenius_norm();
	int rk_svd(double* u, double* sigma, double* v);
};

#endif