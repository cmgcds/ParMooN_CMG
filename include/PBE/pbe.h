#ifndef PBE
#define PBE

#include "hmatrix.h"
#include "sparsematrix.h"

typedef enum {
	kt_rkmatrix = 1,
	kt_convolution = 2,
	kt_lshape = 4,
	kt_hmatrix = 8
} KernelType;

typedef struct _rkdata rkdata;
typedef rkdata* prkdata;

struct _rkdata {
	int r;
	function_1D** fx;
	function_1D** fy;
};

typedef struct _integraloperator integraloperator;
typedef integraloperator* pintegraloperator;

struct _integraloperator {
	int nx;
	double* xgrid;

	KernelType kt;
	prkmatrix rk;
	phmatrix hmatr;
	
	int nr;
	psparsematrix s;
};

pintegraloperator new_integraloperator(int nx, int nr);

void del_integraloperator(pintegraloperator m);

int load_integraloperator(pintegraloperator mo, double* internal_grid, KernelType kt, prkdata lowr, function_2D kl, function_2D kh);

int apply_integraloperator(pintegraloperator mo, double* input, double* output);

int load_hmatrix(pintegraloperator mo, function_2D k);

void generate_mass_matrix(psparsematrix M);

#endif

