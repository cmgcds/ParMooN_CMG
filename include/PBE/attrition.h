#ifndef ATTRITION
#define ATTRITION

#include "sparsematrix.h"
//#include "projection.h"

typedef struct _attrition attrition;
typedef attrition* pattrition;

struct _attrition {
	int nx;

	int nr;
	psparsematrix s;

	function_3D3D* velocity;
};

#define CMAX 9

pattrition new_attrition(int nx, int nr, function_3D3D velocity);

void del_attrition(pattrition m);

//int load_attrition(pattrition mo, double* internal_grid, KernelType kt, prkdata lowr, function_2D kl, function_2D kh);

int apply_attrition(pattrition mo, double* input, double* output, double t);

void generate_mass_matrix(psparsematrix M);

double get_x(int ix, int nx);

double* get_r(int ir, int nr);

int is_on_boundary(double* r);

#endif

