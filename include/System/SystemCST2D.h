/** ************************************************************************ 
*
* @class     TSystemCST2D
* @brief     stores the information of a 2D CST system matrix 
* @author    Jagannath Venkatesan, 
* @date      18.08.15
* @History   Added methods (Jagan, 23.08.14)
 ************************************************************************  */


#ifndef __SYSTEMCST2D__
#define __SYSTEMCST2D__

#include <SquareMatrix2D.h>

/**class for 2D  CST system matrix */
class TSystemCST2D
{
  protected:

    /** DOFs of stress space */    
    int N_S, N_Active, N_DirichletDof;
    
      /** fespace */
    TFESpace2D *FeSpace[5];
      
        /** Fe functions */
    TFEFunction2D *FeFct[3];    

    // used for LPS streamline term
    TFEFunction2D *u1, *u2;
    
    /** Discretization type */
    int Disctype;
    
    /** Solver type */
    int SOLVER;
    
     /** Tensor type */
    int Tensortype;
    
       /** number of matrices in the system matrix*/
    int N_Matrices;
           
    /** Bilinear coefficient   */
    CoeffFct2D *LinCoeffs[1];    

    /** sqstructureA of the system matrix */
    TSquareStructure2D *sqstructure;
  
    /** S is the stiffness/system mat */
    TSquareMatrix2D *SqmatrixS11, *SqmatrixS12, *SqmatrixS21, *SqmatrixS22, *SqmatrixS23, *SqmatrixS32, *SqmatrixS33, *SQMATRICES[8];
    TSquareMatrix **sqmatrices;
  
    TAuxParam2D *NSEaux, *NSEaux_error;
     
    /** Boundary conditon */
    BoundCondFunct2D *BoundaryConditions[3];
  
     /** Boundary values */   
    BoundValueFunct2D *BoundaryValues[3];
        
    /** Discrete form for the equation */
    TDiscreteForm2D *DiscreteFormARhs;

    
  public:
    /** constructor */
     TSystemCST2D(TFESpace2D *stress_fespace, TFEVectFunct2D *Stress, int tensortype, int disctype, int solver, TFESpace2D *Velocity_FeSpace, TFESpace2D* Pressure_FeSpace, TFESpace2D* Deformation_FeSpace = NULL, TFEVectFunct2D *Velocity = NULL);

    /** destructor */
    ~TSystemCST2D();

    /** methods */
    /** Initilize the discrete forms and the matrices */    
    void Init(CoeffFct2D *lincoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *S1BoundValue, BoundValueFunct2D *S2BoundValue, BoundValueFunct2D *S3BoundValue, TAuxParam2D *aux, TAuxParam2D *auxerror);
 
    /** assemble the system matrix */
    void Assemble(double *sol, double *rhs);

    /** solve the system matrix */
    void  Solve(double *sol, double *rhs);
    
    /** measure the error */
    void MeasureErrors(DoubleFunct2D *ExactS1, DoubleFunct2D *ExactS2, DoubleFunct2D *ExactS3, double *s_error);

    
};

#endif
