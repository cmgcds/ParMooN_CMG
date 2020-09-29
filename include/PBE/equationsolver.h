#ifndef MPIMIS_EQUATIONSOLVER
#define MPIMIS_EQUATIONSOLVER

#include "complex.h"

int solveQuadraticEquation(idouble* res, double p_dA, double p_dB, double p_dC);
int solveQuadraticEquationReal(double*  res, int& cRealRoots, double p_dA, double p_dB, double p_dC);
int solveQuadraticEquationReal1(double*  res, int& cRealRoots, double p_dB, double p_dC);
int solveCubicEquation(idouble* res, double p_dB, double p_dC, double p_dD);
int solveCubicEquationReal(double* res, int& cRealRoots, double p_dA, double p_dB, double p_dC, double p_dD);
int solveCubicEquationReal1(double* res, int& cRealRoots, double p_dB, double p_dC, double p_dD);
int solveQuarticEquation(idouble* res, double p_dA, double p_dB, double p_dC, double p_dD, double p_dE);
int solveQuarticEquationReal(double* res, int& cRealRoots, double p_dA, double p_dB, double p_dC, double p_dD, double p_dE);

#endif // MPIMIS_EQUATIONSOLVER
