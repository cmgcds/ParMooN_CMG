// =======================================================================
// @(#)RefTrans1D.h        
//
// Class:      TRefTrans1D
//
// Purpose:    reference transformations for 1D geometric objects
//
// Author:     Sashikumaar Ganesan
//
// History:    17.05.2007 start implementation
// 
// =======================================================================

#include <RefTrans1D.h>

RefTrans1D TRefTrans1D::FindRefTrans1D
        (int N_LocalUsedElements, FE1D *LocalUsedElements)
{
  RefTrans1D rf;

    rf = LineAffin;


  return rf;
}

