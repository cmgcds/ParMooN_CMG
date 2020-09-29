// =======================================================================
// @(#)ItMethod.h        1.6 10/18/99
// 
// Class:       TMultiGridIte
// Purpose:     defines the fixed point iteration
//
// Author:      Volker John
//
// History:     start of implementation 24.10.2000
//
// =======================================================================
#ifndef __MULTIGRIDITE__
#define __MULTIGRIDITE__

#include <ItMethod.h>
#include <NSE_MultiGrid.h>

/** iteration method */
class TMultiGridIte : public TItMethod
{
  protected:
    /** NSE multigrid */  
    TNSE_MultiGrid *mg;
    
    /** start solution */
    int Zero_Start;

  public:
    /** constructor */
    TMultiGridIte(MatVecProc *MatVec, DefectProc *Defect, TItMethod *Prec,
                   int n_aux, int N_Unknowns, TNSE_MultiGrid *MG, int zero_start);

    /** destructor */
    ~TMultiGridIte();
    
    /** iterate routine */
    int Iterate(TSquareMatrix **A, TMatrix **B, double *sol, 
                double *rhs);    
};
#endif
