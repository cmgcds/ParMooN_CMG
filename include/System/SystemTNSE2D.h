/** ************************************************************************ 
*
* @class     TSystemTNSE2D
* @brief     stores the information of a 2D TNSE system matrix 
* @author    Sashikumaar Ganesan, 
* @date      03.09.14
* @History   Added methods (Sashi, 03.09.14)
 ************************************************************************  */

#ifndef __SYSTEMTNSE2D__
#define __SYSTEMTNSE2D__

#include <SystemNSE2D.h>
#ifdef _SMPI
#include <SeqParDirectSolver.h>
#endif

/**class for 2D  TNSE system matrix */
class TSystemTNSE2D : public TSystemNSE2D
{
  protected:
    
    /** M - mass/system mat for TNSE velocity component   */
    TSquareMatrix2D *SqmatrixM11, *SqmatrixM12, *SqmatrixM21, *SqmatrixM22;
  
     /** working rhs, used in AssembleSystMat() */
    double *B;   
   
    /** to store defect */
    double *defect;   
    
    /** factor that multplied with Mat A in working rhs */
    double gamma;      

    /** Discrete form of the M and rhs matrics */
    TDiscreteForm2D *DiscreteFormRhs; 
    
    /** NSE_Rhsaux is used to for assembling rhs only*/
    TAuxParam2D *NSE_Rhsaux;
    
    /** Systmat assemble indicator */
    bool SystMatAssembled;
    
    /** needed for error calculation in time */
    double olderror_l_2_l_2u;
    
#ifdef _SMPI
    TSeqParDirectSolver *P_DS;
#endif     
    
  public:
    /** constructor */
     TSystemTNSE2D(TFESpace2D *velocity_fespace, TFESpace2D *presssure_fespace, TFEVectFunct2D *Velocity, 
                      TFEFunction2D *p, double *sol, double *rhs, int disctype, int nsetype, int solver
#ifdef __PRIVATE__  
                     ,TFESpace2D *Projection_space, TFESpace2D *Stress_FeSpace,   TFESpace2D *Deformation_FeSpace
#endif    
                      );

    /** destrcutor */
    ~TSystemTNSE2D();

    /** methods */
    /** Initilize the discrete forms and the matrices */    
    void Init(CoeffFct2D *lincoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *U1BoundValue, BoundValueFunct2D *U2BoundValue,
              TAuxParam2D *aux, TAuxParam2D *nseaux_error);
    
    
    /** return the stiffness matrix */
    
    
    /** assemble the M, A and rhs */
    void Assemble(double *sol, double *rhs);
        
    /** assemble only the rhs of NSE system */
    void AssembleRhs(double *sol, double *rhs); 
    
    /** scale B matices and assemble rhs based on the \theta scheme  */
    void AssembleSystMat(double scale, double *oldrhs, double *rhs, double *sol);

    /** scale B matices and assemble rhs based on the \theta scheme  */
    void AssembleSystMatNonLinear();
     
    /** restoring the mass matrix */
    void RestoreMassMat();
        
    /** assemble the nonlinear part of the NSE system */
    void AssembleANonLinear(double *sol, double *rhs);

    /** solve the system matrix */
    void  Solve(double *sol);
    
    /** get the resudual of the NSE system */
    void GetTNSEResidual(double *sol, double *res);   
 
    /** measure the error in the NSE */
    void MeasureTNSEErrors(DoubleFunct2D *ExactU1, DoubleFunct2D *ExactU2, DoubleFunct2D *ExactP, double *AllErrors);
};

#endif
