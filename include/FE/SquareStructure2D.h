// =======================================================================
// @(#)SquareStructure2D.h        1.6 09/14/99
//
// Class:       TSquareStructure2D
//
// Purpose:     build and store a structure for a square matrix in 2d
//
// Author:      Gunar Matthies
//
// History:     06.08.1998 start implementation
//
// =======================================================================

#ifndef __SQUARESTRUCTURE2D__
#define __SQUARESTRUCTURE2D__

#include <FESpace2D.h>
#include <SquareStructure.h>

class TSquareStructure2D : public TSquareStructure
{
  protected:
    /** FE space */
    TFESpace2D *FESpace;

  public:
    /** dummy constructor, needed only for derived classes */
    TSquareStructure2D();

    /** generate the matrix structure, only one space needed */
    TSquareStructure2D(TFESpace2D *space);

    /** generate the matrix structure, all arrays are already defined */
    TSquareStructure2D(int n, int N_entries, int *col_ptr,
      int *row_ptr);
    
    /** Generates an empty n*n Structure for a Zero-Matrix */
    explicit TSquareStructure2D(int n);

    /** destructor: free all used arrays */
    ~TSquareStructure2D();

    /** return FESpace */
    TFESpace2D *GetFESpace()
      { return FESpace; }
    
    /** @brief find out if two TSquareStructure2Ds are the same */
    friend bool operator==(const TSquareStructure2D &lhs,
                           const TSquareStructure2D &rhs);
    /** @brief find out if two TSquareStructure2Ds are different */
    friend bool operator!=(const TSquareStructure2D &lhs,
                           const TSquareStructure2D &rhs);
};
#endif
