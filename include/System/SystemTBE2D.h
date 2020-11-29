/** ************************************************************************ 
*
* @class     TSystemTBE2D
* @brief     stores the information of 2D Burgers system matrix 
* @author    Sashikumaar Ganesan, 
* @date      28.11.20
* @History    
 ************************************************************************  */


#ifndef __SYSTEMTBE2D__
#define __SYSTEMTBE2D__

#include <SquareMatrix2D.h>
#include <AssembleMat2D.h>
#include <DirectSparseLinearSolver.h>

/**class for 2D Burgers system matrix */
class TSystemTBE2D
{
  protected:
    
    /** fespace */
    TFESpace2D *FeSpace, *fesp[2];

    /**velo function */
    TFEVectFunct2D *VelocityFct;
    
    /** Fe functions of NSE */
    TFEFunction2D *FeFct[2];    

    /** number of dof in FeSpace */
    int N_DOF, N_Active, N_DirichletDof;
    
    /** Discretization type */
    int Disctype;
    
    /** Solver type */
    int SOLVER;

    /** sparse direct solver */
    TDirectSparseLinearSolver  *DirectSolver;

    /** number of matrices in the system matrix*/
    int N_Matrices;

    /** sqstructure of the system matrix */
    TSquareStructure2D *sqstructure;

    // Assembling */
    TAssembleMat2D *AMatRhsAssemble, *AMatAssembleNonLinear;

    /** M - mass/system mat for TBE velocity component   */
    TSquareMatrix2D *SqmatrixM, *SqmatrixA, *SQMATRICES[2];
    TSquareMatrix **sqmatrices;

    /** working rhs, used in AssembleSystMat() */
    double *RHSs[2], gamma, *B;   
   
    /** to store defect */
    double *defect, olderror_l_2_l_2u; 

    /** Boundary conditon */
    BoundCondFunct2D *BoundaryConditions[2];

     /** Boundary value */   
    BoundValueFunct2D *BoundaryValues[2];
        
    /** Bilinear coefficient   */
    CoeffFct2D *LinCoeffs[1];  

    /** Discrete form for the equation */
    TDiscreteForm2D *DiscreteFormMARhs;
    
    /** Discrete form of A matrics */
    TDiscreteForm2D *DiscreteFormNL; 

    /** BE_Rhsaux is used to for assembling rhs only*/
    TAuxParam2D *BEaux, *BEaux_error;
    
    /** Systmat assemble indicator */
    bool SystMatAssembled;
    
  public:
    /** constructor */
     TSystemTBE2D(TFESpace2D *velocity_fespace,TFEVectFunct2D *Velocity, double *sol, double *rhs, int disctype, int solver);

    /** destrcutor */
    ~TSystemTBE2D();
    
    // /** Initilize the discrete forms and the matrices */
    void Init(CoeffFct2D *lincoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *U1BoundValue,  
                        BoundValueFunct2D *U2BoundValue);
 
    /** assemble the system matrix */
    void Assemble(double *sol, double *rhs);

    /** assemble the system matrix */
    void AssembleA();

    /** assemble the system matrix */
    void AssembleSystMat(double *oldrhs, double *rhs, double *sol);

    void GetTBEResidual(double *sol, double *res);

    /** solve the system matrix */
    void  Solve(double *sol);
    
    /** restore mass mat */
    void RestoreMassMat();

    /** aeesble A mat */
    void AssembleANonLinear(double *sol, double *rhs);

    /** assemble system mat in nonlinear iteration*/
    void AssembleSystMatNonLinear();
    
    /** measure errors */
    void MeasureErrors(DoubleFunct2D *ExactU1, DoubleFunct2D *ExactU2, double *AllErrors);

    // // Return Function for Square Matrix
    // TSquareMatrix2D* ReturnSquareMatrixPointer()
    // {
    //   return sqmatrixA;
    // }
};

#endif
