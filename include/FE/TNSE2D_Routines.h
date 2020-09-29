// ======================================================================
// @(#)TNSE2D_Routines.h        1.5 11/24/99
//
// common declaration for all time dependent Navier-Stokes problems
// ======================================================================

// ======================================================================
// compute turbulent viscosity for LES 
// ======================================================================
double TurbulentViscosity(double hK, double* gradU, double* u, double* uConv);
