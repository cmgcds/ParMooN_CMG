// =======================================================================
// @(#)SquareMatrixNSE2D.h        1.3 11/20/98
// 
// Class:       TSquareMatrixNSE2D
//
// Purpose:     store a square matrix (ansatz = test space) in 2d
//
// Author:      Gunar Matthies
//
// History:     17.11.03 start implementation
//
// =======================================================================

#ifndef __SQUAREMATRIXNSE2D__
#define __SQUAREMATRIXNSE2D__

#include <SquareMatrix2D.h>
#include <SquareStructure2D.h>

class TSquareMatrixNSE2D : public TSquareMatrix2D
{
  protected:
    /** begin array for jb array */
    int *BeginJb;

    /** local number of special edge dof */
    int *jb;

    /** number of dof on each joint */
    int N_DOFperJoint;

    /** coefficient array due to special edge dof */
    double *Alpha;

    /** begin array for bubble correction */
    int *BeginC;

    /** bubble correction coefficients */
    double *C;

    /** begin array for nonconstant pressure */
    int *BeginP;

    /** nonconstant pressure parts */
    double *P;

  public:
    /** generate the matrix */
    TSquareMatrixNSE2D(TSquareStructure2D *squarestructure);

    /** destructor: free Entries array */
    ~TSquareMatrixNSE2D();

    int *GetBeginJb()
    { return BeginJb; }

    int *GetJb()
    { return jb; }

    int GetN_DOFperJoint()
    { return N_DOFperJoint; }

    double *GetAlpha()
    { return Alpha; }

    int *GetBeginC()
    { return BeginC; }

    void SetBeginC(int *beginc)
    { BeginC = beginc; }

    double *GetC()
    { return C; }

    void SetC(double *c)
    { C = c; }

    int *GetBeginP()
    { return BeginP; }

    void SetBeginP(int *beginp)
    { BeginP = beginp; }

    double *GetP()
    { return P; }

    void SetP(double *p)
    { P = p; }
};

#endif
