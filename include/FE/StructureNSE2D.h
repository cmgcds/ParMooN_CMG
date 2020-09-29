// =======================================================================
// @(#)StructureNSE2D.h        1.6 09/14/99
// 
// Class:       TStructureStokes2D
//
// Purpose:     build and store matrix for (Navier-)Stokes equations in 2d
//
// Author:      Gunar Matthies
//
// History:     10.10.2003 start implementation
//
// =======================================================================

#ifndef __STRUCTURENSE2D__
#define __STRUCTURENSE2D__

#include <SquareStructure2D.h>

class TStructureNSE2D : public TSquareStructure2D
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

    void GenerateAlpha();

  public:
    /** generate the matrix structure, only one space needed */
    TStructureNSE2D(TFESpace2D *space);

    /** destructor: free all used arrays */
    ~TStructureNSE2D();

    int *GetBeginJb()
    { return BeginJb; }

    int *GetJb()
    { return jb; }

    int GetN_DOFperJoint()
    { return N_DOFperJoint; }

    double *GetAlpha()
    { return Alpha; }

};

#endif
