// =======================================================================
// @(#)ItMethod.h        1.6 10/18/99
// 
// Class:       TSSORIte
// Purpose:     defines the fixed point iteration
//
// Author:      Volker John
//
// History:     start of implementation 24.10.2000
//
// =======================================================================
#ifndef __SSORITE__
#define __SSORITE__

#include <ItMethod.h>

/** iteration method */
class TSSORIte : public TItMethod
{
  public:
    /** constructor */
  TSSORIte(MatVecProc *MatVec, DefectProc *Defect, TItMethod *Prec,
             int n_aux, int N_Unknowns, int scalar);
  
  /** destructor */
  virtual ~TSSORIte();
  
  /** iterate routine */
  int Iterate(TSquareMatrix **A, TMatrix **B, double *sol, 
              double *rhs);    
};
#endif
