// ======================================================================
// @(#)TNSE3D_Routines.h        1.5 11/24/99
//
// common declaration for all time dependent Navier-Stokes problems
// ======================================================================

// ======================================================================
// compute turbulent viscosity for LES 
// ======================================================================
double TurbulentViscosity3D(double hK, double* gradU, double* u, 
			    double* uConv, double* x, double* y, double* z,
                            double proj_space);

// ======================================================================
// compute stabilization for div--div term
// ======================================================================
double DivDivStab3D(double u1, double u2, double u3, double hK, double eps);  
