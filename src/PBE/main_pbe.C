#include "breakage.h"
#include "aggregation.h"
#include <math.h>
#include <stdlib.h>

int main()
{
	int nxInit = 32, L = 5, nr1 = 17, nr2 = 17, nr3 = 17, nx = (L + 1) * nxInit / 2.0 + 1;
	double Lmin = L_min, Lmax = Lmin * pow(nxInit * pow(2, L - 1) + 1, 1.0/3.0);
	double l1 = 0.01, l2 = 0.01, l3 = 0.5;

	// TO BE CHANGED WITH REAL VELOCITY / SHEAR FLOW
	double* v = allocate_vector(nr1 * nr2 * nr3 * 3);
	fill_vector(v, nr1 * nr2 * nr3 * 3, fabs((double)rand() / (double)20000));
	double* grad_v = allocate_vector(nr1 * nr2 * nr3 * 9);
	fill_vector(grad_v, nr1 * nr2 * nr3 * 9, fabs((double)rand() / (double)10000));

	// TO BE CHANGED WITH REAL INPUT / OUTPUT
	double* input = allocate_vector(nr1 * nr2 * nr3 * nx);
	fill_vector(input, nr1 * nr2 * nr3 * nx, (double)rand() / (double)1000);
	double* output = allocate_vector(nr1 * nr2 * nr3 * nx);
	clear_vector(output, nr1 * nr2 * nr3 * nx);

	pbreakage pb = new breakage(Lmin, Lmax, nxInit, L, nr1, nr2, nr3, l1, l2, l3);
	pb->apply_breakage(input, output, (void*)v);

	paggregation pa = new aggregation(Lmin, Lmax, nxInit, L, nr1, nr2, nr3, l1, l2, l3);
	pa->apply_aggregation(input, output, (void*)grad_v);

	return 0;
}