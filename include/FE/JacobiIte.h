// =======================================================================
// @(#)ItMethod.h        1.6 10/18/99
// 
// Class:       TJacobiIte
// Purpose:     defines the fixed point iteration
//
// Author:      Volker John
//
// History:     start of implementation 24.10.2000
//
// =======================================================================
 #ifndef __JACOBIITE__
 #define __JACOBIITE__

#include <ItMethod.h>
#ifdef _MPI   
   #ifdef __3D__
    #include <ParFECommunicator3D.h>
   #else
    #include <ParFECommunicator2D.h>
   #endif
#endif 

/** iteration method */
class TJacobiIte : public TItMethod
{
  public:
      double *oldSol;
 #ifdef _MPI   
    TParFECommunicator3D *ParComm;
 #endif  
  public:
    /** constructor */
  TJacobiIte(MatVecProc *MatVec, DefectProc *Defect, TItMethod *Prec,
             int n_aux, int N_Unknowns, int scalar
#ifdef _MPI   
                               ,TParFECommunicator3D *ParComm
#endif    
    
  );
  
  /** destructor */
  virtual ~TJacobiIte();
  
  /** iterate routine */
  int Iterate(TSquareMatrix **A, TMatrix **B, double *sol, 
              double *rhs);    
  void Iterate_p(TSquareMatrix **A, TMatrix **B, double *sol, double *rhs
#ifdef _MPI   
                               ,TParFECommunicator3D *ParComm
#endif
                                                   );   
   /** calculate defect */
  void Defect(TSquareMatrix **A, double *sol, double *f, double *d, double &res);
};
 #endif
