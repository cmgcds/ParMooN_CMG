// =======================================================================
// @(#)RefNoRef.C        1.1 10/30/98
//
// Class:       TRefNoRef
// Purpose:     no refinement - only path to shape descriptor
//
// Author:      Volker Behns  31.07.97
//
// =======================================================================

#include <RefNoRef.h>

static Refinements DatEdgeType[MAXN_ORIGEDGES];

// Constructor
TRefNoRef::TRefNoRef(TShapeDesc *shape) : TRefDesc(shape)
{
  int i;

  for (i=0;i<MAXN_ORIGEDGES;i++)
    DatEdgeType[i] = NoRef;

  EdgeType = DatEdgeType;
}

// Methods
