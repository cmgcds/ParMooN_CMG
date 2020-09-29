// =======================================================================
// @(#)SquareMatrixNSE3D.C        1.2 11/20/98
// 
// Class:       TSquareMatrixNSE3D
//
// Purpose:     store a square matrix (ansatz = test space) in 3D
//
// Author:      Gunar Matthies
//
// History:     29.07.05 start implementation
//
// =======================================================================

#include <SquareMatrixNSE3D.h>
#include <StructureNSE3D.h>

TSquareMatrixNSE3D::TSquareMatrixNSE3D(TSquareStructure3D *squarestructure)
  : TSquareMatrix3D(squarestructure)
{
  TStructureNSE3D *structureNSE;

  structureNSE = (TStructureNSE3D *)squarestructure;

  BeginJb = structureNSE->GetBeginJb();
  jb = structureNSE->GetJb();
  N_DOFperJoint = structureNSE->GetN_DOFperJoint();
  Alpha = structureNSE->GetAlpha();

  BeginC = NULL;
  C = NULL;

  BeginP = NULL;
  P = NULL;
}

TSquareMatrixNSE3D::~TSquareMatrixNSE3D()
{
  if(BeginC) delete BeginC;
  if(C) delete C;

  if(BeginP) delete BeginP;
  if(P) delete P;
}
