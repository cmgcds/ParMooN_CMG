// =======================================================================
// @(#)SquareMatrix3D.C        1.2 11/20/98
// 
// Class:       TSquareMatrix3D
//
// Purpose:     store a square matrix (ansatz = test space) in 3d
//
// Author:      Gunar Matthies
//
// History:     10.08.1998 start implementation
//
// =======================================================================

#include <SquareMatrix3D.h>
#include <string.h>

TSquareMatrix3D::TSquareMatrix3D(TSquareStructure3D *squarestructure)
  : TSquareMatrix(squarestructure), structure(squarestructure)
{
}

// not clear! what is the initial value of structure while calling contructor/
// commented by Sashi 
// TSquareMatrix3D::TSquareMatrix3D(int n) 
//  : structure(new TSquareStructure3D(n)), TSquareMatrix(structure)
// {
// }

TSquareMatrix3D::~TSquareMatrix3D()
{
}
