// =======================================================================
// @(#)ItMethod.h        1.6 10/18/99
// 
// Class:       TFixedPointIte
// Purpose:     defines the fixed point iteration
//
// Author:      Volker John
//
// History:     start of implementation 24.10.2000
//
// =======================================================================
#ifndef __FIXEDPOINT__
#define __FIXEDPOINT__

#include <ItMethod.h>
#include <NSE_MultiGrid.h>
#ifdef _MPI   
   #ifdef __3D__
    #include <ParFECommunicator3D.h>
   #else
    #include <ParFECommunicator2D.h>
   #endif
#endif 

/** iteration method */
class TFixedPointIte : public TItMethod
{
  protected: 
#ifdef _MPI   
   #ifdef __3D__
      TParFECommunicator3D *ParComm;
   #else
      TParFECommunicator2D *ParComm;
   #endif
#endif   
  
  public:
   
    /** constructor */
    TFixedPointIte(MatVecProc *MatVec, DefectProc *Defect, TItMethod *Prec,
                   int n_aux, int N_Unknowns, int scalar);

#ifdef _MPI   
    TFixedPointIte(MatVecProc *MatVec, 
                               DefectProc *Defect, 
                               TItMethod *Prec,
                               int n_aux, int n_dof,
                               int scalar, 
  #ifdef  __3D__			       
			       TParFECommunicator3D *parcomm
  #else			       
                               TParFECommunicator2D *parcomm
  #endif
                                );
#endif     
    
    /** destructor */
    virtual ~TFixedPointIte();
    
    /** iterate routine */
    int Iterate(TSquareMatrix **A, TMatrix **B, double *sol, 
                double *rhs);    
};
#endif
