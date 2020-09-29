// =======================================================================
// @(#)Matrix3D.h        1.2 11/20/98
// 
// Class:       TMatrix3D
//
// Purpose:     store a  matrix (ansatz != test space)
//
// Author:      Gunar Matthies
//
// History:     26.08.1998 start implementation
//
// =======================================================================

#ifndef __MATRIX3D__
#define __MATRIX3D__

#include <Structure3D.h>
#include <Matrix.h>

class TMatrix3D : public TMatrix
{
  protected:
    /** matrix structure */
    TStructure3D *structure;

  public:
    /** generate the matrix */
    TMatrix3D(TStructure3D *structure);

    /** destructor: free Entries array */
    ~TMatrix3D();

    TStructure3D *GetStructure()
    { return structure; }
    
    /** @brief set all Dirichlet rows to zero. That means all rows where the 
     * test space has nonactive degrees of freedom. 
     */
    void resetNonActive();
};

#endif
