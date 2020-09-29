 // =======================================================================
// @(#)SquareMatrix2D.h        
// 
// Class:       TSquareMatrix2D
//
// Purpose:     store a square matrix (ansatz = test space) in 1d
//
// Author:      Sashikumaar Ganesan
//
// History:     17.05.2007 start implementation
//
// =======================================================================

#ifndef __SQUAREMATRIX1D__
#define __SQUAREMATRIX1D__

#include <SquareMatrix.h>
#include <SquareStructure1D.h>

class TSquareMatrix1D : public TSquareMatrix
{
  protected:
    /** matrix strcuture */
    TSquareStructure1D *structure;

  public:
    /** generate the matrix */
    TSquareMatrix1D(TSquareStructure1D *squarestructure);

    /** destructor: free Entries array */
    ~TSquareMatrix1D();

    /** return FESpace */
    TFESpace1D *GetFESpace() const
    { return structure->GetFESpace(); }

    /** return used matrix structure */
    TSquareStructure1D *GetMatrixStructure() const
    { return structure; }

};

#endif

