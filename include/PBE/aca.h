#ifndef ACA
#define ACA

typedef enum {
  HLIB_ACA_DEFAULT,		/* Try to bound relative Euclidean error */
  HLIB_ACA_ABSOLUTE,		/* Try to bound absolute Frobenius error */
  HLIB_ACA_RELATIVE,		/* Try to bound relative Frobenius error */
  HLIB_ACA_ABSOLUTE2,		/* Try to bound absolute Euclidean error */
  HLIB_ACA_RELATIVE2		/* Try to bound relative Euclidean error */
} ACAStrategy;

/* Compute an adaptive cross approximation of a matrix with modified pivoting
   (ACA+ algorithm). The return value will be the rank k <= kmax, the first
   k colums of A and B will be filled by the rank-k-approximation of the matrix.
   A: Row matrix, containing at least rows*kmax double variables
   B: Column matrix, containing at least cols*kmax double variables
   rows: Number of rows
   cols: Number of columns
   row_off: Row offset, rows [row_off ... row_off+rows-1] are considered
   col_off: Column offset, columns [col_off ... col_off+cols-1] are considered
   entry: Callback function for computing a matrix entry
   data: Additional data for callback function "entry"
   strategy: Truncation strategy
   kmax: Maximal rank
   eps: Tolerance for truncation */
int
newaca_fill_block(double *A, double *B, int rows, int cols,
	       int row_off, int col_off,
	       double (*entry)(int row, int col, void *data), void *data,
	       int kmax, double eps, ACAStrategy strategy);

/* Compute an adaptive cross approximation of a matrix. The return value will
   be the rank k <= kmax, the first k colums of A and B will be filled by
   the rank-k-approximation of the matrix.
   A: Row matrix, containing at least rows*kmax double variables
   B: Column matrix, containing at least cols*kmax double variables
   rows: Number of rows
   cols: Number of columns
   fillrow: Callback function for computing a matrix row
   fillcol: Callback function for computing a matrix column
   data: Additional data for callback functions "fillrow" and "fillcol"
   strategy: Truncation strategy
   kmax: Maximal rank
   eps: Tolerance for truncation */
int
aca_fill_block(double *A, double *B, int rows, int cols,
	       int row_start, int col_start,
	       double (*entry)(int row, int col, void *data),
	       void *data, int kmax, double eps, ACAStrategy strategy);

#endif

