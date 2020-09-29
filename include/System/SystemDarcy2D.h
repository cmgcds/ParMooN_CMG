/** ************************************************************************ 
*
* @class     TSystemDarcy2D
* @brief     stores the information of a 2D Darcy system matrix 
* @author    Ulrich Wilbrandt
* @date      15.03.15
 ************************************************************************  */


#ifndef __SYSTEMDARCY2D__
#define __SYSTEMDARCY2D__

#include <SquareMatrix2D.h>
#include <LocalAssembling2D.h>

/**class for 2D scalar system matrix */
class TSystemDarcy2D
{
  protected:
    
    /** fespaces for velocity and pressure */
    TFESpace2D *fe_spaces[2];
    
    /** number of matrices in the system matrix */
    int N_Matrices;

    // ( A  B1' )   ( 0 2 )
    // ( B2 C   )   ( 3 1 )
    
    /** A is the stiffness/system mat for stationary problem   */
    TSquareMatrix2D *sq_matrices[2];
    
    TMatrix2D *rect_matrices[2];
    
    /** Boundary conditon (one for u.n and one for pressure) */
    BoundCondFunct2D *BoundaryConditions[2];

     /** Boundary value */ 
    BoundValueFunct2D *BoundaryValues[2];
    
  public:
    /** constructor */
     TSystemDarcy2D(TFESpace2D **fespaces);

    /** destrcutor */
    ~TSystemDarcy2D();
    
    /** Initilize the discrete forms and the matrices */
    void Init(BoundCondFunct2D **BoundCond, BoundValueFunct2D **BoundValue);
 
    /** assemble the system matrix */
    void Assemble(LocalAssembling2D& la, double *sol, double *rhs);

    /** solve the system matrix */
    void  Solve(double *sol, double *rhs);
};

#endif // __SYSTEMMATDARCY2D__
