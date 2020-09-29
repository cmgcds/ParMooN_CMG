// =======================================================================
// @(#)BoundPart.C        1.1 10/30/98
//
// Class:       TBoundPart
// Purpose:     a closed part of boundary
//
// Author:      Volker Behns  28.06.97
//
// =======================================================================

#include <BoundPart.h>
#include <MooNMD_Io.h>

// Constuctor
TBoundPart::TBoundPart (int n_comps)
{
  N_BdComps = n_comps;
#ifdef __2D__
  BdComps = new TBoundComp2D*[N_BdComps];
#else
  BdComps = new TBoundComp3D*[N_BdComps];
#endif
}

// Methods
