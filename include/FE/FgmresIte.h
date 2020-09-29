// =======================================================================
// @(#)ItMethod.h        1.6 10/18/99
// 
// Class:       TFgmresIte
// Purpose:     defines the fixed point iteration
//
// Author:      Volker John
//
// History:     start of implementation 24.10.2000
//
// =======================================================================
#ifndef __FGMRESITE__
#define __FGMRESITE__

#ifdef _MPI   
   #ifdef __3D__
    #include <ParFECommunicator3D.h>
   #else
    #include <ParFECommunicator2D.h>
   #endif
#endif 
      
#include <ItMethod.h>

/** iteration method */
class TFgmresIte : public TItMethod
{
  protected:

  /** arrays for flexible gmres depending on restart */
  double *s;
  double *cosi;
  double *sn;

  /** arrays for flexible gmres depending on number of dof */
  double *d;
  double *z;
  
  /** matrices for flexible gmres depending on restart */
  double **H;
  double **v;
  double **zv;
  
  int scalar;

#ifdef _MPI   
   #ifdef __3D__
      TParFECommunicator3D *ParComm,*ParCommU,*ParCommP;
   #else
      TParFECommunicator2D *ParComm,TParFECommunicator2D *parcommP;
   #endif
#endif  
  
  public:
    /** constructor */
  //  TFgmresIte(MatVecProc *MatVec, DefectProc *Defect, TItMethod *Prec,
    //           int n_aux, int N_Unknowns, int scalar);
 
//#ifdef _MPI
    TFgmresIte(MatVecProc *MatVec, DefectProc *Defect, TItMethod *Prec,
               int n_aux, int N_Unknowns, int scalar
#ifdef _MPI
  #ifdef  __3D__			       
			       ,TParFECommunicator3D *parcommU=NULL,
			TParFECommunicator3D *parcommP = NULL
  #else			       
                               ,TParFECommunicator2D *parcommU=NULL,
			TParFECommunicator2D *parcommP = NULL
  #endif
#endif
	      );
    
//     TFgmresIte(MatVecProc *MatVec,
// 			DefectProc *Defect,
// 			TItMethod *Prec,
// 			int n_aux, int n_dof,
// 			int scalar_,
//   	#ifdef  __3D__
// 			       TParFECommunicator3D *parcomm
// #else		       
// 				TParFECommunicator2D *parcomm
// 	#endif 
//                       );
//     
// #endif
    /** destructor */
    virtual ~TFgmresIte();
    
    /** iterate routine */
    int Iterate(TSquareMatrix **A, TMatrix **B, double *sol, 
                double *rhs);    
    
    double Dotprod(int x, double* v1, double* v2);
    
    void update(double* &v1);
};
#endif
