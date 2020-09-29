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
#ifndef __MULTIGRIDSCAITE__
#define __MULTIGRIDSCAITE__

#include <ItMethod.h>

#ifdef __2D__
#include <MultiGrid2D.h>
#else
#include <MultiGrid3D.h>
#endif

/** iteration method */
class TMultiGridScaIte : public TItMethod
{
  protected:
    /** NSE multigrid */  
#ifdef __2D__
    TMultiGrid2D *mg;
#else
    TMultiGrid3D *mg;
#endif

    /** start solution */
    int Zero_Start;

  public:
    /** constructor */
#ifdef __2D__
    TMultiGridScaIte(MatVecProc *MatVec, DefectProc *Defect, TItMethod *Prec,
                     int n_aux, int N_Unknowns, TMultiGrid2D *MG,
                     int zero_start);
#else
    TMultiGridScaIte(MatVecProc *MatVec, DefectProc *Defect, TItMethod *Prec,
                     int n_aux, int N_Unknowns, TMultiGrid3D *MG,
                     int zero_start);
#endif

    /** destructor */
    ~TMultiGridScaIte();
    
    /** iterate routine */
    int Iterate(TSquareMatrix **A, TMatrix **B, double *sol, 
                double *rhs);
};
#endif
