/** ************************************************************************ 
* @brief     source file for TDirectSparseLinearSolver
* @author    Sashikumaar Ganesan
* @date      28.11.20
* @History 
 ************************************************************************  */

#include <DirectSparseLinearSolver.h>
#include <MainUtilities.h>
#include <DirectSolver.h>
#include <Database.h>
#include "stdlib.h"
#include <LinAlg.h>

#ifdef _MPI
  #include <mkl.h>
  #include <omp.h>
#endif


#include <fstream>

extern "C"
{
#include "umfpack.h"
}

TDirectSparseLinearSolver::TDirectSparseLinearSolver(TSquareMatrix *matrix, int n_rhs)
{
  Matrix = matrix;
  N_Rhs = n_rhs;
  // Rhs = rhs;
  N_Eqn = Matrix->GetN_Columns();

  int *Row_orig = Matrix->GetRowPtr();
  int *KCol_orig = Matrix->GetKCol();
  double *Values_orig = Matrix->GetEntries();

  N_Entries = Row_orig[N_Eqn];
  KCol = new int[N_Entries];
  Row = new int[N_Eqn + 1];
  Values = new double[N_Entries];

  memcpy(Values, Values_orig, N_Entries * SizeOfDouble);
  memcpy(KCol, KCol_orig, N_Entries * SizeOfInt);
  memcpy(Row, Row_orig, (N_Eqn + 1) * SizeOfInt);  

  // check ordering of the matrix
  if (matrix->GetColOrder() != 1)
   {
    // sort matrix
    OutPut("umfpack: reordering of the columns needs to be performed" << endl);
    exit(0);
   }

  double t1 = GetTime();
  ret = umfpack_di_symbolic(N_Eqn, N_Eqn, Row, KCol, Values, &Symbolic, NULL, NULL);
  double t2 = GetTime();

  // error occured
  if (ret != 0)
  {
    OutPut("error in umfpack_di_symbolic " << ret << endl);
    //exit(4711);
  }
  
}

TDirectSparseLinearSolver::~TDirectSparseLinearSolver()
{
 
  umfpack_di_free_symbolic(&Symbolic);

  delete[] Values;
  delete[] KCol;
  delete[] Row;  
}

/** solver_flag = 0 => factorize, solve and delete */
/** solver_flag = 1 => factorize and solve  */
/** solver_flag = 2 => solve  */

void TDirectSparseLinearSolver::DirectSolve(int solver_flag, double *rhs, double *sol)
{
  double *Values_orig = Matrix->GetEntries();
  memcpy(Values, Values_orig, N_Entries * SizeOfDouble);

  if (solver_flag == 0 || solver_flag == 1)
  {
    ret = umfpack_di_numeric(Row, KCol, Values, Symbolic, &Numeric, NULL, NULL);
    // umfpack_di_free_symbolic(&Symbolic);
    //t3 = GetTime();
  
    // error occured
    if (ret != 0)
    {
    OutPut("error in umfpack_di_numeric " << ret << endl);
    exit(4711);
    }
  }

  for (int i = 0; i < N_Rhs; i++)
  {
   //Sol = sol[i];
   Sol = sol + i * N_Eqn;
   Rhs = rhs + i * N_Eqn;

   ret = umfpack_di_solve(UMFPACK_At, Row, KCol, Values, Sol, Rhs, Numeric, NULL, NULL);

    //t4 = GetTime();
    if (ret != 0)
    {
      OutPut("error in umfpack_di_solve " << ret << endl);
      //exit(4711);
    }
  } // for(i=0; i<N_Rhs; i++)

  if (solver_flag == 0)
  {
    umfpack_di_free_numeric(&Numeric);
  }

//  cout <<"DirectSolve" <<endl;
//  exit(0);

} // void TDirectSparseLinearSolver::DirectSolve()
