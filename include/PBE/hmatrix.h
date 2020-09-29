#ifndef HMATRIX
#define HMATRIX

#include "rkmatrix.h"
#include "fullmatrix.h"

typedef enum {
	hs_strong, 
	hs_weak
} HStrategy;

typedef struct _hmatrix hmatrix;
typedef hmatrix *phmatrix;

struct _hmatrix {
	int rows;
	int cols;
	int block_rows;
	int block_cols;
	int l;
	prkmatrix r;
	pfullmatrix f;
	phmatrix* sons;
};

phmatrix new_hmatrix(int r, int c);
void del_hmatrix(phmatrix h);

int create_hmatrix(phmatrix h, double (*entry)(int row, int col, void* data), void* data, int row_start, int col_start);
int create_rkmatrix(phmatrix h, double (*entry)(int row, int col, void* data), void* data, int row_start, int col_start);
int create_fullmatrix(phmatrix h, double (*entry)(int row, int col, void* data), void* data, int row_start, int col_start);

int hmvm(phmatrix h, double* b, double* c);
int hmmm(phmatrix h, double* B, int B_cols, double* C);
int trihmvm(phmatrix h, char* uplo, double* b,double* c);
int trihmmm(phmatrix h, char* uplo, double* B,int B_cols, double* C);

#endif

