// =======================================================================
// @(#)InnerEdge.h  
// 
// Class:       TInnerEdge
// Purpose:     class for inner edges in 3D
//
// Author:      Sashikumaar Ganesan  03.09.2010
//
// History:
//
// ======================================================================= 

#ifndef __INNEREDGE__
#define __INNEREDGE__

#include <Edge.h>

/** an edge in a 3D grid */
class TInnerEdge : public TEdge
{
  protected:

  public:

  /** constructor with neighbours */
   TInnerEdge(int n_Neibs, TBaseCell **neighbs);

};

#endif
