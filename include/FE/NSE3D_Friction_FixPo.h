// ======================================================================
// @(#)NSE2D_FixPo.h        1.2 04/13/00
//
// common declaration for all Navier-Stokes problems
// ======================================================================

#ifndef __NSE3D_FRICTION_FIXPO__
#define __NSE3D_FRICTION_FIXPO__

#include <Enumerations.h>

// ======================================================================
// Type 3, Standard Galerkin, with friction, D(u):D(v)
// ======================================================================
void NSType3GalerkinDD3DFriction(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Galerkin, friction, D(u):D(v)
// ======================================================================
void NSType3_4NLGalerkinDD3DFriction(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

#endif
