// =======================================================================
// @(#)MumpsSolver.h
//
// Class:      TMumpsSolver
// Purpose:    Solve equation system by  MPI based MUMPS parallel direct solvers
//
// Author:     Sashikumaar Ganesan (30.09.09)
//
// History:    Start of implementation 30.09.09 (Sashikumaar Ganesan)
//
// =======================================================================
#if defined(_MPI) || defined(_SMPI)
 extern "C"
 {
  #  include "dmumps_c.h"
 }

#ifndef __MUMPSSOLVER__
#define __MUMPSSOLVER__

 #include <SquareMatrix.h>
 #include <Matrix.h>
 #include <SquareMatrix3D.h>
 #include <Matrix3D.h>
 #include <ParFECommunicator3D.h>

/** general class for all parallel direct solvers */

class TMumpsSolver
{
  protected:

  /** MPI_Comm for which the fespace communications are needed */
  MPI_Comm Comm;

  /**  internal pointers of MUMPS solver*/
  DMUMPS_STRUC_C id;
 
  /** global rhs */
  double *MumpsRhs;
 
  /** */
  bool FactorFlag;
  
  
  public:
   /** constructor */
   TMumpsSolver(int N_Eqns, int M_dist_Nz, int *M_dist_Irn, int *M_dist_Jcn, int N_Rhs);

   void FactorizeAndSolve(double *Mat_loc, double *rhs);

   void Solve(double *Mat_loc, double *rhs);
   
   void Clean();
  
    /** destructor */
    ~TMumpsSolver();
};
#endif


#endif
