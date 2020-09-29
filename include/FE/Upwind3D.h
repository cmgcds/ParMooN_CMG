// =======================================================================
// @(#)Upwind.h        1.4 10/18/99
//
// Purpose:     upwind stabilization
//              do upwind assembling for first order nonconforming elements
//
// Authors:     Volker Behns / Gunar Matthies  18.10.99
//
// =======================================================================


#ifndef __UPWIND3D__
#define __UPWIND3D__

void UpwindForNavierStokes3D(TSquareMatrix3D *sqmatrix, TFEFunction3D *u1,
                             TFEFunction3D *u2, TFEFunction3D *u3);

void UpwindForConvDiff(TSquareMatrix3D *sqmatrix, double *RHS,
                       TFESpace3D *fespace, TDiscreteForm3D
                       *DiscreteForm);

#endif
