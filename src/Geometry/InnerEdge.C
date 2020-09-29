// =======================================================================
// @(#)InnerEdge  
// 
// Class:       TInnerEdge
// Purpose:     class for inner edges in 3D
//
// Author:      Sashikumaar Ganesan  03.09.2010
//
// History:
//
// ======================================================================= 

#include <InnerEdge.h>

/** constructor with neighbours */
TInnerEdge::TInnerEdge(int n_Neibs, TBaseCell **neighbs):TEdge(n_Neibs, neighbs)
{
  EdgeID = InnerEdge;
}


