/** ************************************************************************ 
*
* @class     TSystemTCD3D
* @brief     stores the information of a timedependent part of a 3D scalar system matrix 
* @author    Sashikumaar Ganesan
* @date      24.01.15
* @History 
 ************************************************************************  */


#ifndef __SYSTEMTCD3D__
#define __SYSTEMTCD3D__

#include <SquareMatrix3D.h>
#include <SystemCD3D.h>

/**class for 3D scalar system matrix */
class TSystemTCD3D : public TSystemCD3D
{
  protected:
#ifdef _MPI
    TParDirectSolver *TDS;
#endif   
    /** M mass matrix */
    TSquareMatrix3D **sqmatrixM;
    
    /** working rhs, used in AssembleSystMat() */
    double *B;   
   
    /** to store defect */
    double *defect;   
    
    /** factor that multplied with Mat A in working rhs */
    double gamma;   
    
    bool factorize;

    
    /** instance of the Assemble class */
    TAssembleMat3D **MMatRhsAssemble;
    
//     /** Stiffness part of the SUPG matrix */
//     TSquareMatrix3D *sqmatrixK;    
//     
//     /** time-consistent part of the SUPG matrix */
//     TSquareMatrix3D *sqmatrixS;
    
    /** Discrete form of the M and rhs matrics */
    TDiscreteForm3D *DiscreteFormMRhs, *DiscreteFormRhs; 
    
    /** Systmat assemble indicator */
    bool SystMatAssembled;
    
  public:
    /** constructor */
     TSystemTCD3D(int N_levels, TFESpace3D **fespaces, double **sol, double **rhs, int disctype, int solver);

    /** destrcutor */
    ~TSystemTCD3D();

    /** methods */
    void Init(CoeffFct3D *BilinearCoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *BoundValue, TAuxParam3D *aux);
    
    /** return the stiffness matric */
    TSquareMatrix3D **GetAMatrix()
    { return sqmatrixA; }
    
    /** assemble the Mass mat and rhs */
    void AssembleMRhs(); 
    
    /** assemble the stifness mat and rhs */
    void AssembleARhs();   
    
//     /** M = M + (tau*THETA1)*A */ 
//     /** B = (tau*THETA1)*rhs +(tau*THETA2)*oldrhs + [ M - (tau*THETA2)A]*oldsol */  
    void AssembleSystMat(double *oldrhs, double *oldsol, double *rhs, double *sol
#ifdef _MPI
                         , double **Rhs_array
#endif
                         );
    
    /** restoring the mass matrix */
    void RestoreMassMat();
    
     /** solve the system matrix */
    void Solve(double *sol);  
    
    /** return the residual of the system for the given sol*/
    double GetResidual(double *sol);
    
    double value(double *sol,int N){
      int i;
      double sum=0.;
      for(i=0;i<N;i++)	sum+=sol[i];
      return sum;
    }
    
};

#endif

