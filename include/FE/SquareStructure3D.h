// =======================================================================
// @(#)SquareStructure3D.h        1.6 09/14/99
// 
// Class:       TSquareStructure3D
//
// Purpose:     build and store a structure for a square matrix in 3D
//
// Author:      Gunar Matthies
//
// History:     06.08.1998 start implementation
//
// =======================================================================

#ifndef __SQUARESTRUCTURE3D__
#define __SQUARESTRUCTURE3D__

#include <FESpace3D.h>
#include <SquareStructure.h>

class TSquareStructure3D : public TSquareStructure
{
  protected:
    /** FE space */
    TFESpace3D *FESpace;

  public:
    /** dummy constructor, needed only for derived classes */
    TSquareStructure3D();
  
    /** generate the matrix structure, only one space needed */
    TSquareStructure3D(TFESpace3D *space);

    /** generate the matrix structure, all arrays are already defined */
    TSquareStructure3D(int n, int N_entries, int *col_ptr,
      int *row_ptr);
    
    /** Generates an empty n*n Structure for a Zero-Matrix */
    explicit TSquareStructure3D(int n);

    /** destructor: free all used arrays */
    ~TSquareStructure3D();

    /** return FESpace */
    TFESpace3D *GetFESpace()
    { return FESpace; }

};

#endif
