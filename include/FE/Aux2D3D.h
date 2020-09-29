// =======================================================================
// %W% %G%
//
// Class:       TAux2D3D
// Purpose:     calculate parameters from 3d fefuntion for the use
//              in 2d
//
// Author:      Gunar Matthies (15.05.01)
//
// History:     start of implementation 15.05.01 (Gunar Matthies)
// =======================================================================

#ifndef __AUX2D3D__
#define __AUX2D3D__

#include <FEFunction3D.h>

/** calculate parameters from 3d fefuntion for the use in 2d */
class TAux2D3D
{
  protected:
    /** collection of 3d cells */
    TCollection *Coll;

    /** 3d cell number */
    int *CellNumbers;

    /** joint numbers for 3d cells */
    int *JointNumbers;

    /** 3d finite element function */
    TFEFunction3D *FEFunct;

    /** 3d finite element space */
    TFESpace3D *FESpace;

    /** BeginIndex */
    int *BeginIndex;

    /** GlobalNumbers */
    int *GlobalNumbers;

    /** value vector of FEFunct */
    double *Values;

    /** shift in parameter list */
    int Shift;

  public:
    /** constructor */
    TAux2D3D(int *cellnumbers, int *jointnumbers, TFEFunction3D *fefunct,
             int shift);

    /** calculate gradient for local coordinates (xi,eta) on face
        JointNumbers[num] of cell CellNumbers[num] */
    void GetGradient(int num, int N_Points, double *xi, double *eta,
                     double **Param);
};

#endif
