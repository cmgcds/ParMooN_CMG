/** ************************************************************************ 
*
* @class     TSystemTCD2D
* @brief     stores the information of a timedependent part of a 2D scalar system matrix 
* @author    Sashikumaar Ganesan
* @date      08.08.14
* @History 
 ************************************************************************  */


#ifndef __SYSTEMTCD2D__
#define __SYSTEMTCD2D__

#include <SquareMatrix2D.h>
#include <SystemCD2D.h>
#include <CDSystemTimeDG.h>
#include <CDSystemTimeDG_1.h>

/**class for 2D scalar system matrix */
class TSystemTCD2D : public TSystemCD2D
{
  protected:
        
    /** M mass matrix */
    TSquareMatrix2D *sqmatrixM;
    
    /** working rhs, used in AssembleSystMat() */
    double *B;   
   
    /** to store defect */
    double *defect;   
    
    /** factor that multplied with Mat A in working rhs */
    double gamma;   
    
    /** Stiffness part of the SUPG matrix */
    TSquareMatrix2D *sqmatrixK;    
    
    /** time-consistent part of the SUPG matrix */
    TSquareMatrix2D *sqmatrixS;
    
    /** Discrete form of the M and rhs matrics */
    TDiscreteForm2D *DiscreteFormMRhs; 
    
    /** Systmat assemble indicator */
    bool SystMatAssembled;
    
   /** dG time steppings */
   TCDSystemTimeDG *TimeDG;
//    TCDSystemTimeDG_1 *TimeDG_1;
  
  public:
    /** constructor */
     TSystemTCD2D(TFESpace2D *fespace, int disctype, int solver);

    /** destrcutor */
    ~TSystemTCD2D();

    /** methods */
    void Init(CoeffFct2D *BilinearCoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *BoundValue);
    
    /** return the stiffness matric */
    TSquareMatrix2D *GetAMatrix()
    { return sqmatrixA; }
    
    /** assemble the Mass mat and rhs */
    void AssembleMRhs(TAuxParam2D *aux, double *sol, double *rhs); 
    
    /** assemble the stifness mat and rhs */
    void AssembleARhs(TAuxParam2D *aux, double *sol, double *rhs);   
    
    /** M = M + (tau*TDatabase::TimeDB->THETA1)*A */ 
    /** B = (tau*TDatabase::TimeDB->THETA1)*rhs +(tau*TDatabase::TimeDB->THETA2)*oldrhs + [ M - (tau*TDatabase::TimeDB->THETA2)A]*oldsol */  
    void AssembleSystMat(double *oldrhs, double *oldsol, double *rhs, double *sol);
    
    /** restoring the mass matrix */
    void RestoreMassMat();
    
     /** solve the system matrix */
    void Solve(double *sol, double *rhs);  
    
    /** return the residual of the system for the given sol*/
    double GetResidual(double *sol);
    
};

#endif
