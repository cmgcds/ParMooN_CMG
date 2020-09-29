// =======================================================================
// @(#)RefTriBis0Desc.h        1.1 10/30/98
//
// Class:       TRefTriBis0Desc
// Purpose:     refinement descriptor for bisection of a triangle
//              bisection of edge 0
//
// Author:      Volker Behns  18.07.97
//
// =======================================================================

#ifndef __REFTRIBIS0DESC__
#define __REFTRIBIS0DESC__

#include <RefDesc.h>

#define TRIBI0MAXN_VpC    3
#define TRIBI0MAXN_CpV    2
#define TRIBI0MAXN_EpC    3
#define TRIBI0MAXN_CpE    2
#define TRIBI0MAXN_EpV    3
#define TRIBI0MAXN_iVpE   1
#define TRIBI0MAXN_nVpoE  3
#define TRIBI0MAXN_nEpoE  2
#define TRIBI0N_E         3

/** refinement descriptor for bisection of a triangle */
class TRefTriBis0Desc : public TRefDesc
{
  public:
    // Constructor
    /** build the descriptor for refining a triangle by bisecting edge 0 */
    TRefTriBis0Desc(TShapeDesc *shape);

    // Methods
};

#endif
