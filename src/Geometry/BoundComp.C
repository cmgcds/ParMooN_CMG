// =======================================================================
// @(#)BoundComp.C        1.2 08/03/99
//
// Class:       TBoundComp
// Purpose:     prototype of a component of a boundary part
//
// Author:      Volker Behns  18.06.97
//
// =======================================================================

#include <BoundComp.h>

// Constructor
TBoundComp::TBoundComp(int id, int ref)
{
  ID = id;
  refID = ref;
  FreeBoundaryStatus = false;
}

// Methods
