#include "sparsematrix.h"
#include <malloc.h>

psparsematrix new_sparsematrix(int r, int c, int cmax)
{
	psparsematrix s;
	s = (psparsematrix)malloc(sizeof(sparsematrix));

	s->rows = r;
	s->cols = c;
	s->cmax = cmax;

	s->nzindices = (int*)malloc(cmax * c * sizeof(int));
	s->data = allocate_vector(cmax * c);

	return s;
}

void del_sparsematrix(psparsematrix s)
{
	s->rows = 0;
	s->cols = 0;
	s->cmax = 0;
	free((void*)s->nzindices);
	free_vector(s->data);
	s = 0;
}

void matrix_times_sparsematrix(psparsematrix s, double* A, int A_rows, double* B)
{
	int i, j, h;
	for(i = 0; i < s->cols; i++)
		for(j = 0; j < A_rows; j++)
		{
			B[i * A_rows + j]  = 0;
			for(h = 0; h < s->cmax; h++)
			{
				if(s->nzindices[s->cmax * i + h] >= 0)
					B[i * A_rows + j] += s->data[s->cmax * i + h] * A[s->nzindices[s->cmax * i + h] * A_rows + j];
			}
		}
}

