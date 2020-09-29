#ifndef MPIMIS_BREAKAGE
#define MPIMIS_BREAKAGE

#include "basictools.h"
#include "constants.h"
#include "equationsolver.h"
//#include "newton.h"
#include "timer.h"

typedef class Breakage breakage;
typedef breakage* pbreakage;

class Breakage {
public:
	// Constructor for uniform grid
	Breakage(double _xmin, double _xmax, int _na, int _nx, int _ny, int _nz) {
			 xmin = pow(_xmin, 3), xmax = pow(_xmax, 3), na = _na, nx = _nx, ny = _ny, nz = _nz;
			 a = allocate_vector(na); g = allocate_vector(na); lg = (int*)malloc((na + 1) * sizeof(int)); 
			 freq_integral_coeff = allocate_vector(na); plus_integral_coeff = allocate_vector(na);
			 a = allocate_vector(na); fill_uniform_grid(xmin, xmax, na, a); transform_grid_m2l(xmin, xmax, na, a);
			 fill_freq_integral_coeff(); };

	~Breakage() { free_vector(a);
		free_vector(freq_integral_coeff); free_vector(plus_integral_coeff); free_vector(g); free((void*)lg); };

private:
	double xmin;
	double xmax;
	double* a;	int na;
	int* lg; double* g;
	double* freq_integral_coeff;
	double* plus_integral_coeff;

	int nx, ny, nz;

private:
	// Plus and minus parts with respect to one single spatial point
	// sizoef(input / output) = na, sizeof(data = velocity) = 3
	void breakage_plus(double* input, double* output, double* data);
	void breakage_minus(double* input, double* output);

	int is_on_boundary(int ir, double* vr);

	double get_value_at_x(double* fi, double x);

	void fill_freq_integral_coeff();

	int solvequarticequation(double* roots, int& cRealRoots, double C_r, double cube_x);

//	double lower_bound_g(double& x, double& CE);
	void fill_g(double& CE);
	void fill_lg(double& CE);
	void fill_plus_integral_coeff();
	double definite_integral(double& fa, double& fb, double& a, double& b, double c, double& x1, double& x2);
	double frequency(double* input);

public:
	// Breakage operator with the whole input and output
	// sizoef(input / output) = nx * ny * nz * na
	// sizeof(data = velocity) = 3 * nx * ny * nz
	void apply_breakage(double* input, double* output, double* v);
};

#endif // MPIMIS_BREAKAGE
