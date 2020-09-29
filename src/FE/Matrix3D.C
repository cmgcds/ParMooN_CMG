// =======================================================================
// @(#)Matrix3D.C        1.2 11/20/98
// 
// Class:       TMatrix3D
//
// Purpose:     store a  matrix3D (ansatz != test space)
//
// Author:      Gunar Matthies
//
// History:     26.08.1998 start implementation
//
// =======================================================================

#include <Matrix3D.h>
#include <string.h>

TMatrix3D::TMatrix3D(TStructure3D *structure)
 : TMatrix(structure)
{
  this->structure = structure;
}

TMatrix3D::~TMatrix3D()
{
}


void TMatrix3D::resetNonActive()
{
  int n_active = this->structure->GetTestSpace()->GetN_DegreesOfFreedom()
                -this->structure->GetTestSpace()->GetN_Dirichlet();
  int * rowPtr = this->structure->GetRowPtr();
  int index_nonactive = rowPtr[n_active];
  int n_nonactive_entries = rowPtr[this->structure->GetN_Rows()]
                           - index_nonactive;
  memset(Entries + index_nonactive, 0.0, n_nonactive_entries * SizeOfDouble);
}
