// =======================================================================
// @(#)BDEdge3D.h  
// 
// Class:       TBDEdge3D
// Purpose:     class for boundary edges in 3D
//
// Author:      Sashikumaar Ganesan  03.09.2010
//
// History:
//
// ======================================================================= 

#ifndef __BDEDGE3D__
#define __BDEDGE3D__

#include <Edge.h>


/** an edge in a 3D grid */
class TBDEdge3D : public TEdge
{
  protected:

  public:

    /** constructor with neighbours */
    TBDEdge3D(int n_Neibs, TBaseCell **neighbs);
#ifdef _MPI
  void set_Bd_id(int key)
  {
    Bd_id = key;
  }
  
  int get_Bd_id()
  {      
    return Bd_id;
  }
      
#endif

};

#endif
