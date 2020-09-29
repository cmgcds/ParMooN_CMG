#ifndef SPARSEMATRIX
#define SPARSEMATRIX

#include "basictools.h"

typedef struct _sparsematrix sparsematrix;
typedef sparsematrix *psparsematrix;
typedef const sparsematrix *pcsparsematrix;

struct _sparsematrix {
	int rows;
	int cols;
	int cmax; // maximum possible number of entries for each column

	int* nzindices; // indices of rows of nonzero entries for each column 
					// (cmax * cols dimensional array). If -1 then 
	double* data; // values of nonzero entries (cmax * cols dimensional array)
};

psparsematrix new_sparsematrix(int r, int c, int cmax);
void del_sparsematrix(psparsematrix s);

void matrix_times_sparsematrix(psparsematrix s, double* A, int A_rows, double* B);

#endif

