// =======================================================================
// @(#)SquareMatrix2D.C        1.2 11/20/98
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

#include <SquareMatrix1D.h>
#include <string.h>

TSquareMatrix1D::TSquareMatrix1D(TSquareStructure1D *squarestructure)
  : TSquareMatrix(squarestructure), structure(squarestructure)
{
}

TSquareMatrix1D::~TSquareMatrix1D()
{
}
 
