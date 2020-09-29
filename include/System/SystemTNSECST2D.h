 
/** ************************************************************************ 
*
* @class     TSystemTNSECST2D
* @brief     stores the information of a monolithic 2D T-NSE-CST system matrix 
* @author    Jagannath Venkatesan, 
* @date      15.07.16
* @History   
 ************************************************************************  */


#ifndef __SYSTEMTNSECST2D__
#define __SYSTEMTNSECST2D__

    
#include <SystemNSECST2D.h>

/**class for 2D T-NSE-CST system matrix */
class TSystemTNSECST2D : public TSystemNSECST2D
{
  protected:

   /** MU - mass/system mat for TNSE velocity component   */
    TSquareMatrix2D *SqmatrixMU11, *SqmatrixMU12, *SqmatrixMU21, *SqmatrixMU22;
 
    /** MS - mass/system mat for TCST stress component   */
    TSquareMatrix2D *SqmatrixMS11, *SqmatrixMS12, *SqmatrixMS21, *SqmatrixMS22, *SqmatrixMS23, *SqmatrixMS32, *SqmatrixMS33;
    
    /** working rhs, used in AssembleSystMat() */
    double *B;
    
    /** used to handle newton-type nonlinearity   */ 
    double *NLRhs;
    
    /** to store defect */
    double *defect;  
    
    /** factor that multplied with Mat A in working rhs */
    double gamma;
    
    /** Systmat assemble indicator */
    bool SystMatAssembled;
    
    /** needed for error calculation in time */
    double olderror_l_2_l_2u, olderror_l_2_l_2s, olderror_l_2_h_1u, olderror_l_2_h_1s, olderror_l_2_l_2p, olderror_l_2_h_1p;
    
    
  public:
    /** constructor */
     TSystemTNSECST2D(int N_FESpaces, TFESpace2D **FE_Spaces, int N_FEFunctions, TFEFunction2D **FeFunctions, int disctype_NSECST, int solver);

    /** destructor */
    ~TSystemTNSECST2D();

     /** methods */
    /** Initilize the discrete forms and the matrices */ 
    void Init(CoeffFct2D *lincoeffs, BoundCondFunct2D **BoundCond, BoundValueFunct2D **BoundValue,  
			   TAuxParam2D *TNSECSTaux, TAuxParam2D *TNSEauxerror, TAuxParam2D *TCSTauxerror);
       
    /** assemble the M, system matrices and rhs */
   void Assemble(double *sol, double *rhs);
   
   /** Assemble the non-linear matrices only A,G,D **/ 
   void AssembleNonLinear(double *sol, double *rhs);
   
   /** assemble only the rhs */
    void AssembleRhs(double *sol, double *rhs);
   
   /** scale the matices and assemble rhs based on the \theta scheme  */
    void AssembleSystMat(double scale, double *oldrhs, double *rhs, double *sol);
    
    /** scale the matices and assemble lhs and rhs based on the \theta scheme  */
    void AssembleSystMatNonLinear();
    
    
    /** restoring the mass matrix */
    void RestoreMassMat();
    
    /** solve the system matrix */
    void  Solve(double *sol);
    
    /** get the residual */
    void GetTNSECSTResidual(double *sol, double *res);
    
        /** measure the errors */
    void MeasureTNSECSTErrors(DoubleFunct2D *ExactU1, DoubleFunct2D *ExactU2, DoubleFunct2D *ExactP, DoubleFunct2D *ExactS1, DoubleFunct2D *ExactS2, DoubleFunct2D *ExactS3, double *AllErrors );
};

#endif
