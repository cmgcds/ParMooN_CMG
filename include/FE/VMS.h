// =======================================================================
// VMS.h
//
// Purpose:     routines for projection-based VMS
//
// Author:       Volker John  2006/05/18
//
// =======================================================================

#ifndef __VMS__
#define __VMS__

#include <SquareMatrix2D.h>
#include <Matrix2D.h>

#ifdef __3D__
  #include <SquareMatrix3D.h>
  #include <Matrix3D.h>
#endif

#ifdef __2D__
void VMSProjectionUpdateMatrices(int N_U,int N_Active,int N_L,
				 TSquareMatrix2D **SQMATRICES, 
				 TMatrix2D **MATRICES);

void LumpMassMatrixToDiag(TSquareMatrix2D *M);

#endif // __2D__

#ifdef __3D__
void VMS_ProjectionUpdateMatrices(int N_U,int N_Active,int N_L,
                             TSquareMatrix3D **SQMATRICES, TMatrix3D **MATRICES);

void VMS_ProjectionExplUpdateRhs(int N_U, int N_Active, int N_L, TFEVectFunct3D *u,
				 TSquareMatrix3D **SQMATRICES, TMatrix3D **MATRICES, 
				 double *rhs_vms_expl);

void LumpMassMatrixToDiag(TSquareMatrix3D *M);

void ComputeVMSProjection(TMatrix3D *G11, TMatrix3D *G22, TMatrix3D *G33,
			  TSquareMatrix3D *MatrixL, TFEFunction3D *u_1, 
			  TFEFunction3D *u_2, TFEFunction3D *u_3,
			  TFEVectFunct3D *vms_projection_fe);

void ComputeSizeOfSmallScales(TMatrix3D *matG11, TMatrix3D *matG22, TMatrix3D *matG33,
			  TSquareMatrix3D *MatrixL, TFEFunction3D *u1, 
			  TFEFunction3D *u2, TFEFunction3D *u3,
			  TFEVectFunct3D *vms_projection_fe, double *size_small_scales);

void MeanAndLargestSize( TFESpace3D *projection_space, double *size_small_scales, 
                         double *mean, double *largest_size);

void AdaptProjectionSpace(TFESpace3D *projection_space, 
			  double *size_small_scales, 
			  FE3D  *fes, 
			  double mean, 
			  double mean_time_average, 
			  double largest_size, 
			  double max_time_average, 
			  double *label_space);


#endif // __3D__ 

#endif

