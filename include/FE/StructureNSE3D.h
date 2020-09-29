// =======================================================================
// @(#)StructureNSE3D.h        1.6 09/14/99
// 
// Class:       TStructureStokes3D
//
// Purpose:     build and store matrix for (Navier-)Stokes equations in 3D
//
// Author:      Gunar Matthies
//
// History:     29.07.2005 start implementation
//
// =======================================================================

#ifndef __STRUCTURENSE3D__
#define __STRUCTURENSE3D__

#include <SquareStructure3D.h>

class TStructureNSE3D : public TSquareStructure3D
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
    TStructureNSE3D(TFESpace3D *space);

    /** destructor: free all used arrays */
    ~TStructureNSE3D();

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
