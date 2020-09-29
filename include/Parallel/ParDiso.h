// =======================================================================
// @(#)ParDiso.h
//
// Class:      TParDiso
// Purpose:    Solve equation system by ParDiso routines
//
// Author:     Sashikumaar Ganesan & Abdus Shamim (01.05.15)
//
// History:    Start of implementation 01.05.15 (Sashikumaar Ganesan & Abdus Shamim) 
//
// =======================================================================

#ifdef _OMPONLY

#ifndef __PARDISO__
#define __PARDISO__

#include <Database.h>
// #include <SquareMatrix.h>
// #include <Matrix.h>
// #include <SquareMatrix3D.h>
// #include <Matrix3D.h>
// #include <ParFECommunicator3D.h>

class TParDiso
{
  protected:
    
  void *pt[64];
  int iparam[64];
  double dparam[64];
  int phase, nrhs, Nmax_system, matrix_number, matrix_type, N_Eqns, N_Nz;
  //csr format
  int *RowPtr, *KCol;
  int Solver;
  int perm_user, msglvl, ierror;
  
  double idum, ddum;
  
  double *rhsptr,*solptr;
  
  public:
   /** constructor */
   TParDiso(int neqns, int nnz, int* rowptr, int* kcol);
   
   void FactorizeAndSolve(double *Mat, double *rhs, double *sol,bool Factorize);

   void Solve(double *Mat, double *rhs, double *sol);
   
   void Clean(double *Mat);
  
    /** destructor */
    ~TParDiso();
};
#endif


#endif
