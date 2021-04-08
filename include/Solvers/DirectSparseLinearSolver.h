/** ************************************************************************ 
*
* @class     TDirectSparseLinearSolver
* @brief     base class for Direct Sparse Linear Solver
* @author    Sashikumaar Ganesan 
* @date      28.11.20
* @History    
 ************************************************************************  */


#ifndef __DIRECTSPARSELINEARSOLVER__
#define __DIRECTSPARSELINEARSOLVER__

#include <AllClasses.h>
#include <Constants.h>

#ifdef _MPI

#endif

#ifdef _OMPONLY

#endif

/** class for Direct Sparse Linear Solver */
class TDirectSparseLinearSolver
{
  protected:  

  /** UMFPACK variables */
  void *Symbolic, *Numeric;
  int ret;

  /** system info */
  double *Sol, *Rhs, *Values;

  /** matrix info */
  int N_Eqn, N_Rhs, N_Entries, N_Rhs_Disp, *KCol, *Row;

  /** system matrix */
  TSquareMatrix *Matrix;


  public:
    /** Constructor*/
    TDirectSparseLinearSolver(TSquareMatrix *matrix, int n_rhs);
    
    /** destrcutor */
    ~TDirectSparseLinearSolver();
    
   /** methods */
   void DirectSolve(int solver_flag, double *rhs, double *sol);


};

#endif
