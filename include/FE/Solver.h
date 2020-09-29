// =======================================================================
// @(#)Solver.h        1.9 06/27/00
// 
// Purpose:     solve equation system
//
// Author:      Gunar Matthies (17.08.98)
//
// History:     start of implementation 17.08.98 (Gunar Matthies)
//
// =======================================================================

#ifndef __SOLVER__
#define __SOLVER__

#include <SquareMatrix.h>
#include <Matrix.h>

/** solve equation system */

void Solver(TSquareMatrix *matrix, double *rhs, double *sol);

void Solver(TSquareMatrix *matrix, double *rhs, double *sol, int N_Rhs);

void Solver(int what, TSquareMatrix *matrix, double *rhs, double *sol,
            int N_Rhs);

void Solver(TSquareMatrix *matrix, TMatrix *matrixmortar,
            double *rhs, double *sol);

void Solver(TSquareMatrix *matrixa11, TMatrix *matrixa12,
            TMatrix *matrixa21, TSquareMatrix *matrixa22,
            double *rhs1, double *rhs2,
            double *sol1, double *sol2);

void Solver(TSquareMatrix *matrixa11, TSquareMatrix *matrixa12,
            TSquareMatrix *matrixa22, TSquareMatrix *matrixa23,
            TSquareMatrix *matrixa32, TSquareMatrix *matrixa33,
            double *rhs1, double *rhs2, double *rhs3,
            double *sol1, double *sol2, double *sol3);

/*******************************************************************/
/* CONNECT SYSTEM FOR VMM [KL02]                                   */
/*******************************************************************/
void Solver(TSquareMatrix *matrixa11, TMatrix *matrixa12, TMatrix *matrixa13,
            TMatrix *matrixa21, TMatrix *matrixa31,  TSquareMatrix *matrixa22,
            double *rhs1, double *rhs2, double *rhs3,
            double *sol1, double *sol2, double *sol3);

/*******************************************************************/
/* Rosenbrock                                                      */
/*******************************************************************/
void Solver(TSquareMatrix *matrixa11, TSquareMatrix *matrixa12,
            TSquareMatrix *matrixa21, TSquareMatrix *matrixa22,
            double *rhs1, double *rhs2,
            double *sol1, double *sol2);

void Solver(TSquareMatrix *matrixa11, TSquareMatrix *matrixa12,
            TSquareMatrix *matrixa22, TSquareMatrix *matrixa23,
            TSquareMatrix *matrixa32, TSquareMatrix *matrixa33,
            double *rhs1, double *rhs2, double *rhs3,
            double *sol1, double *sol2, double *sol3);
            
void Solver(TSquareMatrix *matrixa11, TSquareMatrix *matrixa12, TSquareMatrix *matrixa13,
            TSquareMatrix *matrixa21, TSquareMatrix *matrixa22, TSquareMatrix *matrixa23,
            TSquareMatrix *matrixa31, TSquareMatrix *matrixa32, TSquareMatrix *matrixa33,
            double *rhs1, double *rhs2, double *rhs3,
            double *sol1, double *sol2, double *sol3);

void Solver(TSquareMatrix **sqmatrices, TMatrix **matrices, double *rhs,
            double *sol, MatVecProc *MatVect, DefectProc *Defect,
            TMultiGrid2D *MG, int N_Unknowns, int ns_type);

#ifdef __2D__
/** STOKES TYPE SADDLE POINT PROBLEM (NSTYPE 1) */
void Solver(TSquareMatrix *matrixA, TMatrix *matrixB1,
            TMatrix *matrixB2, double *rhs, double *sol);

/** STOKES TYPE SADDLE POINT PROBLEM (NSTYPE 2) */
void Solver(TSquareMatrix *matrixA, TMatrix *matrixB1,
            TMatrix *matrixB2, TMatrix *matrixB3,
            TMatrix *matrixB4, double *rhs, double *sol);

/** STOKES TYPE SADDLE POINT PROBLEM (NSTYPE 3) */
void Solver(TSquareMatrix *matrixA11,TSquareMatrix *matrixA12,
            TSquareMatrix *matrixA21,TSquareMatrix *matrixA22,
            TMatrix *matrixB1, TMatrix *matrixB2,
            double *rhs, double *sol);

/** STOKES TYPE SADDLE POINT PROBLEM (NSTYPE 4) */
void Solver(TSquareMatrix *matrixA11,TSquareMatrix *matrixA12,
            TSquareMatrix *matrixA21,TSquareMatrix *matrixA22,
            TMatrix *matrixB1, TMatrix *matrixB2,
            TMatrix *matrixB3, TMatrix *matrixB4, 
            double *rhs, double *sol);

/** BRAESS--SARAZIN SMOOTHER (NSTYPE 1) */
void Solver(TSquareMatrix *matrixA, TMatrix *matrixB1,
            TMatrix *matrixB2, double *rhs, double *sol,
            int para0);

/** BRAESS--SARAZIN SMOOTHER (NSTYPE 2) */
void Solver(TSquareMatrix *matrixA, TMatrix *matrixB1,
            TMatrix *matrixB2, TMatrix *matrixB3,
            TMatrix *matrixB4, double *rhs, double *sol,
            int para0);

/** BRAESS--SARAZIN SMOOTHER (NSTYPE 3) */
void Solver(TSquareMatrix *matrixA11,TSquareMatrix *matrixA12,
            TSquareMatrix *matrixA21,TSquareMatrix *matrixA22,
            TMatrix *matrixB1, TMatrix *matrixB2,
            double *rhs, double *sol, int para0);

/** BRAESS--SARAZIN SMOOTHER (NSTYPE 4) */
void Solver(TSquareMatrix *matrixA11,TSquareMatrix *matrixA12,
            TSquareMatrix *matrixA21,TSquareMatrix *matrixA22,
            TMatrix *matrixB1, TMatrix *matrixB2,
            TMatrix *matrixB3, TMatrix *matrixB4, 
            double *rhs, double *sol, int para0);

/**  AUXILIARY PROBLEM FOR VASSILEVSKI/LAZAROV APPROACH */
void Solver(TSquareMatrix *matrixA, TMatrix *matrixB1,
            TMatrix *matrixB2, double *rhs, double *sol,double delta);

/** STOKES TYPE SADDLE POINT PROBLEM WITH MORTAR (NSTYPE 2) */
void Solver(TSquareMatrix *matrixA, TMatrix *matrixB1,
            TMatrix *matrixB2, TMatrix *matrixB3,
            TMatrix *matrixB4, TMatrix *matrixmortar,
            double *rhs, double *sol);
#endif
#ifdef __3D__
/** STOKES TYPE SADDLE POINT PROBLEM (NSTYPE 1) */
void Solver(TSquareMatrix *matrixA, TMatrix *matrixB1,
            TMatrix *matrixB2, TMatrix *matrixB3, 
            double *rhs, double *sol);

/** BRAESS--SARAZIN SMOOTHER (NSTYPE 1)     */
void Solver(TSquareMatrix *matrixA, TMatrix *matrixB1,
            TMatrix *matrixB2, TMatrix *matrixB3, 
            double *rhs, double *sol,
            int para0);

#endif
#endif
