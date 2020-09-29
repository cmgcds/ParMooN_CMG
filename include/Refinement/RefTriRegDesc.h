// =======================================================================
// @(#)RefTriRegDesc.h        1.1 10/30/98
//
// Class:       TRefTriRegDesc
// Purpose:     refinement descriptor for regular refinement of a triangle
//
// Author:      Volker Behns  17.07.97
//
// =======================================================================

#ifndef __REFTRIREGDESC__
#define __REFTRIREGDESC__

#include <RefDesc.h>

#define TRIRRN_E         3
#define TRIRRMAXN_VpC    3
#define TRIRRMAXN_CpV    3
#define TRIRRMAXN_EpC    3
#define TRIRRMAXN_CpE    2
#define TRIRRMAXN_EpV    4
#define TRIRRMAXN_iVpE   1
#define TRIRRMAXN_nVpoE  3
#define TRIRRMAXN_nEpoE  2

/** refinement descriptor for regular refinement of a triangle */
class TRefTriRegDesc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a triangle */
    TRefTriRegDesc(TShapeDesc *shape);

    // Methods
};

#endif
