#include "breakage.h"
#include "aggregation.h"
#include "BrAgg.h"
#include <math.h>
#include <stdlib.h>

int main()
{
	int n = 9;
	int nx = 1, ny = 1, nz = 1, na = ((n - 1) / 2.0) * (cLevel + 1) + 1;

	srand(time(NULL));

	double* v;
	v = allocate_vector(nx * ny * nz * 3);
	fill_vector(v, nx * ny * nz * 3, 1);

	double* grad_v;
	grad_v = allocate_vector(nx * ny * nz * 9);
	fill_vector(grad_v, nx * ny * nz * 9, 1);

	double* input; input = allocate_vector(nx * ny * nz * na);
	fill_vector(input, nx * ny * nz * na, 1);

	double* output; output = allocate_vector(nx * ny * nz * na);
	clear_vector(output, nx * ny * nz * na);

	double* temp; temp = allocate_vector(nx * ny * nz);
	fill_vector(temp, nx * ny * nz, 1);

	double L_max = 1e-3, f_max = 1e3;
	Breakage_agglomeration(nx, ny, nz, n, input, v, grad_v, temp, output, L_max, f_max);

	return 0;
}
