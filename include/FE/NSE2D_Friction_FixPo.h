// ======================================================================
// @(#)NSE2D_FixPo.h        1.2 04/13/00
//
// common declaration for all Navier-Stokes problems
// ======================================================================

#ifndef __NSE2D_FRICTION_FIXPO__
#define __NSE2D_FRICTION_FIXPO__

#include <Enumerations.h>

// ======================================================================
// Type 3, Standard Galerkin with friction, D(u):D(v)
// ======================================================================
void NSType3GalerkinDDFriction(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Galerkin with local friction, D(u):D(v)
// ======================================================================
void NSType3GalerkinDDFrictionLocal(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, with friction, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void NSType3UpwindFriction(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, SDFEM, with friction, (grad u, grad v)
// ======================================================================
void NSType4SDFEMFriction(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, SDFEM, with friction, (grad u, grad v)
// ======================================================================
void NSType4SDFEMFrictionRST(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, SDFEM, with friction, D(u):D(v)
// ======================================================================
void NSType4SDFEMDDFriction(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, SDFEM, with friction on rhs, D(u):D(v)
// ======================================================================
void NSType4SDFEMDDFrictionRhs(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Galerkin with friction, D(u):D(v)
// ======================================================================
void NSType3_4NLGalerkinDDFriction(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Standard Galerkin with local friction, D(u):D(v)
// ======================================================================
void NSType3_4NLGalerkinDDFrictionLocal(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void NSType3_4NLUpwindFriction(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, SDFEM, with friction, (grad u, grad v)
// ======================================================================
void NSType4NLSDFEMFriction(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, SDFEM, with friction, (grad u, grad v), Lutz
// ======================================================================
void NSType4NLSDFEMFrictionRST(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, SDFEM, with friction, D(u):D(v)
// ======================================================================
void NSType4NLSDFEMDDFriction(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, SDFEM, with friction on rhs, D(u):D(v)
// ======================================================================
void NSType4NLSDFEMDDFrictionRhs(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

#endif
