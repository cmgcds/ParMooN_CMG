// =======================================================================
// @(#)SquareMatrixNSE2D.C        1.2 11/20/98
// 
// Class:       TSquareMatrixNSE2D
//
// Purpose:     store a square matrix (ansatz = test space) in 2d
//
// Author:      Gunar Matthies
//
// History:     17.11.03 start implementation
//
// =======================================================================

#include <SquareMatrixNSE2D.h>
#include <StructureNSE2D.h>

TSquareMatrixNSE2D::TSquareMatrixNSE2D(TSquareStructure2D *squarestructure)
  : TSquareMatrix2D(squarestructure)
{
  TStructureNSE2D *structureNSE;

  structureNSE = (TStructureNSE2D *)squarestructure;

  BeginJb = structureNSE->GetBeginJb();
  jb = structureNSE->GetJb();
  N_DOFperJoint = structureNSE->GetN_DOFperJoint();
  Alpha = structureNSE->GetAlpha();

  BeginC = NULL;
  C = NULL;

  BeginP = NULL;
  P = NULL;
}

TSquareMatrixNSE2D::~TSquareMatrixNSE2D()
{
  if(BeginC) delete BeginC;
  if(C) delete C;

  if(BeginP) delete BeginP;
  if(P) delete P;
}
