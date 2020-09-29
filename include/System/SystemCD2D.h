/** ************************************************************************ 
*
* @class     TSystemCD2D
* @brief     stores the information of a 2D scalar system matrix 
* @author    Sashikumaar Ganesan, 
* @date      08.08.14
* @History   Added methods (Sashi, 22.08.14)
 ************************************************************************  */


#ifndef __SYSTEMCD2D__
#define __SYSTEMCD2D__

#include <SquareMatrix2D.h>

/**class for 2D scalar system matrix */
class TSystemCD2D
{
  protected:
    
    /** fespace */
    TFESpace2D *FeSpace;

    /** number of dof in FeSpace */
    int N_DOF, N_Active;
    
    /** Discretization type */
    int Disctype;
    
    /** Solver type */
    int SOLVER;
       
    /** number of matrices in the system matrix*/
    int N_Matrices;

    /** sqstructure of the system matrix */
    TSquareStructure2D *sqstructure;

    /** A is the stiffness/system mat for stationary problem   */
    TSquareMatrix2D *sqmatrixA, *SQMATRICES[3];;
    TSquareMatrix **sqmatrices;
    
    /** Boundary conditon */
    BoundCondFunct2D *BoundaryConditions[1];

     /** Boundary value */   
    BoundValueFunct2D *BoundaryValues[1];
        
    /** Discrete form for the equation */
    TDiscreteForm2D *DiscreteFormARhs;
 
    
  public:
    /** constructor */
     TSystemCD2D(TFESpace2D *fespace, int disctype, int solver);

    /** destrcutor */
    ~TSystemCD2D();
    
    /** Initilize the discrete forms and the matrices */
    void Init(CoeffFct2D *BilinearCoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *BoundValue);
 
    /** assemble the system matrix */
    void Assemble(TAuxParam2D *aux, double *sol, double *rhs);

    /** solve the system matrix */
    void  Solve(double *sol, double *rhs);
    
    /** mass and volume */
    void GetMassAndArea(TFEFunction2D *fefunction, double *parameters);  

    // Return Function for Square Matrix
    TSquareMatrix2D* ReturnSquareMatrixPointer()
    {
      return sqmatrixA;
    }
};

#endif
