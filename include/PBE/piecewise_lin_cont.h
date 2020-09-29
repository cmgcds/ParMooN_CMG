#include "basictools.h"
#include "fepc_easy.h"
#include "faltung.h"

class piecewise_constant
{
public:
	int n;
	double* grid; // array of size n + 1 // grid[0] = a, grid[n] = b
	double* value; // array of size n

	double integral();
};

class piecewise_linear
{
public:
	int n;
	double* grid; // array of size n + 1 // grid[0] = a, grid[n] = b
	double* start_value; // array of size n
	double* end_value; // array of size n

	double integral();
};

class piecewise_linear_continuous
{
public:
	piecewise_linear_continuous() {};

	int n;
	double* grid; // array of size n + 1 // grid[0] = a, grid[n] = b
	double* coeff; // array of size n + 1 // coefficiets of the basis functions (hat functions). The same as the values at gridpoints.

	double L1();
	double integral();
	double mass();

	piecewise_linear_continuous& operator = (const piecewise_linear_continuous& f);

	int multiplication_with_function (double* f, piecewise_constant& w);
	int multiplication_with_function (double* f, piecewise_linear& w);

	int ProjectToMassConservedExact();
	int ProjectToMassConservedShift();
	int ProjectToMassConservedPSD();

	int IsContinuousProjectionOf(piecewise_constant& f, double* diag, double* subdiag);
	int IsContinuousProjectionOf(piecewise_linear& f, double* diag, double* subdiag);

	piecewise_constant DiscontinuosRepresentation();
};

int fill_diagonal_part_of_Gramm_matrix(int n, double* grid, double* diag);
int fill_subdiagonal_part_of_Gramm_matrix(int n, double* grid, double* subdiag);
