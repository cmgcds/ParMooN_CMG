#ifndef IFUNCTIONS
#define IFUNCTIONS

#include "constants.h"
#include "basictools.h"

double nu_frag(double x, double* r, function_3D3D velocity, double t);
double V_attr(double x, double* r, function_3D3D velocity, double t);
double E_kin(double x, double* r, function_3D3D velocity, double t);
const double L_min = 32 * mu * Gamma_Kr / (3 * H_V * H_V);
double L_max(double x, double* r, function_3D3D velocity, double t);
double L_star(double x, double* r, function_3D3D velocity, double t);
double b(double* r, function_3D3D velocity, double t);

#endif

