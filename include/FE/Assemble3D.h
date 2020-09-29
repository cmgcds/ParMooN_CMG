// =======================================================================
// %W% %G%
// 
// Purpose:     assemble matrix and right-hand side
//
// Author:      Gunar Matthies (10.08.98)
//
// History:     start of implementation 10.08.98 (Gunar Matthies)
//
// =======================================================================

#ifndef __ASSEMBLE3D__
#define __ASSEMBLE3D__

#include <AllClasses.h>
#include <Constants.h>

/** a function from a finite element space */
void Assemble3D(int n_fespaces, TFESpace3D **fespaces,
                int n_sqmatrices, TSquareMatrix3D **sqmatrices,
                int n_matrices, TMatrix3D **matrices,
                int n_rhs, double **rhs, TFESpace3D **ferhs,
                TDiscreteForm3D *DiscreteForm,
                BoundCondFunct3D **BoundaryConditions,
                BoundValueFunct3D **BoundaryValues,
                TAuxParam3D *parameters);


/** a function from a finite element space */
void Assemble3DSlipBC(int n_fespaces, TFESpace3D **fespaces,
                int n_sqmatrices, TSquareMatrix3D **sqmatrices,
                int n_matrices, TMatrix3D **matrices,
                int n_rhs, double **rhs, TFESpace3D **ferhs,
                TDiscreteForm3D *DiscreteForm,
                BoundCondFunct3D **BoundaryConditions,
                BoundValueFunct3D **BoundaryValues,
                TAuxParam3D *parameters);

void ModifyMatrixSlipBC(TSquareMatrix3D **sqmatrices, TMatrix3D **matrices,
			int N_U, double *rhs);

/** assemble mixed finite elements such as Raviart-Thomas or
 * Brezzi-Douglas-Marini.
 */
void Assemble3D_mixed(int n_fespaces, TFESpace3D **fespaces,
int n_sqmatrices, TSquareMatrix3D **sqmatrices,
int n_matrices, TMatrix3D **matrices,
int n_rhs, double **rhs, TFESpace3D **ferhs,
TDiscreteForm3D *DiscreteForm3D,
BoundCondFunct3D **BoundaryConditions,
BoundValueFunct3D **BoundaryValues,
TAuxParam3D *Parameters);



#endif // __ASSEMBLE3D__
