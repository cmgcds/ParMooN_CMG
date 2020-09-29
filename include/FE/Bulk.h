// =======================================================================
// Buld.h
//
// Purpose:   common routines for bulk precipitation in 2d/3d and 3d/4d
//
// Author:    Volker John  
//
// =======================================================================

#ifndef __BULK__
#define __BULK__

// ======================================================================
// lump matrix to diagonal matrix
// the sparsity pattern of the matrix is not condensed
// ======================================================================

#ifdef __2D__
void LumpMassMatrixToDiag_Bulk(TSquareMatrix2D *M);
#endif    
#ifdef __3D__
void LumpMassMatrixToDiag_Bulk(TSquareMatrix3D *M);
#endif   

double calculate_dp_50(int N, double *size, double *number);

#endif

