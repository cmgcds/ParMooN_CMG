// =======================================================================
// @(#)RefTrans2D.C        1.1 10/30/98
//
// Class:      TRefTrans2D
//
// Purpose:    reference transformations for 2D geometric objects
//
// Author:     Gunar Matthies
//
// History:    08.07.97 start implementation
// 
// =======================================================================

#include <RefTrans2D.h>

RefTrans2D TRefTrans2D::FindRefTrans2D
        (int N_LocalUsedElements, FE2D *LocalUsedElements)
{
  RefTrans2D rf;

  if( ((int)(LocalUsedElements[0])) < 5)
    rf = TriaAffin;
  else
    if( ((int)(LocalUsedElements[0])) < 10)
      rf = QuadAffin;
    else
      rf = QuadBilinear;

  return rf;
}
