// =======================================================================
// @(#)Upwind.h        1.4 10/18/99
//
// Purpose:     upwind stabilization
//              do upwind assembling for first order nonconforming elements
//
// Authors:     Volker Behns / Gunar Matthies  18.10.99
//
// =======================================================================


#ifndef __UPWIND__
#define __UPWIND__

void UpwindForNavierStokes(CoeffFct2D *Coeff, TSquareMatrix2D *sqmatrix, 
			   TFEFunction2D *u1, TFEFunction2D *u2);

void UpwindForConvDiff(CoeffFct2D *Coeff, 
		       TSquareMatrix2D *sqmatrix, double *RHS,
                       TFESpace2D *fespace, TDiscreteForm2D *DiscreteForm,
		       TFEFunction2D *u1, TFEFunction2D *u2,
		       int ConvIsVelo);

/******************************************************************************/
//
// IMPROVED MIZUKAMI-HUGHES METHOD (Knobloch, CMAME 2007)
//
/******************************************************************************/

void ComputeParametersMizukamiHughes(TBaseCell *cell, int cell_no, 
				     TFEFunction2D *u, CoeffFct2D *Coeffs,
				     BoundCondFunct2D *BoundaryCondition,
                                     int *dof, int ActiveBound,
				     double *c_mh);

void MizukamiHughes(TSquareMatrix2D *sqmatrix, 
		    double *RHS,
		    TFESpace2D *fespace, 
		    TFEFunction2D *u, 
		    CoeffFct2D *Coeffs,
		    BoundCondFunct2D *BoundaryCondition);


#endif
