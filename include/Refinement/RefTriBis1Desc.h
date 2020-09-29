// =======================================================================
// @(#)RefTriBis1Desc.h        1.1 10/30/98
//
// Class:       TRefTriBis1Desc
// Purpose:     refinement descriptor for bisection of a triangle
//              bisection of edge 1
//
// Author:      Volker Behns  18.07.97
//
// =======================================================================

#ifndef __REFTRIBIS1DESC__
#define __REFTRIBIS1DESC__

#include <RefDesc.h>

#define TRIBI1MAXN_VpC    3
#define TRIBI1MAXN_CpV    2
#define TRIBI1MAXN_EpC    3
#define TRIBI1MAXN_CpE    2
#define TRIBI1MAXN_EpV    3
#define TRIBI1MAXN_iVpE   1
#define TRIBI1MAXN_nVpoE  3
#define TRIBI1MAXN_nEpoE  2
#define TRIBI1N_E         3

/** refinement descriptor for bisection of a triangle */
class TRefTriBis1Desc : public TRefDesc
{
  public:
    // Constructor
    /** build the descriptor for refining a triangle by bisecting edge 1 */
    TRefTriBis1Desc(TShapeDesc *shape);

    // Methods
};

#endif
