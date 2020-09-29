// =======================================================================
// @(#)BDEdge3D 
// 
// Class:       TBDEdge3D
// Purpose:     class for boundary edges in 3D
//
// Author:      Sashikumaar Ganesan  03.09.2010
//
// History:
//
// ======================================================================= 

#include <BDEdge3D.h>

/** constructor with neighbours */
TBDEdge3D::TBDEdge3D(int n_Neibs, TBaseCell **neighbs):TEdge(n_Neibs, neighbs)
{
  EdgeID = BDEdge3D;
}



