#ifndef FULLMATRIX
#define FULLMATRIX

#include "basic.h"

typedef struct _fullmatrix fullmatrix;
typedef fullmatrix *pfullmatrix;

struct _fullmatrix {
	int rows;
	int cols;
	double* data;	
};

pfullmatrix new_fullmatrix(int r, int c);
void del_fullmatrix(pfullmatrix f);

void fullmatrix_times_matrix(pfullmatrix f, EvalMode mode, double* B, int n, double* C, double alpha, double beta);
void fullmatrix_times_vector(pfullmatrix f, EvalMode mode, double* b, double* c, double alpha, double beta);

#endif

