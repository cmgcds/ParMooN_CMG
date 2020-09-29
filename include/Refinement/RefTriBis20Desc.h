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

#ifndef __REFTRIBIS20DESC__
#define __REFTRIBIS20DESC__

#include <RefDesc.h>

#define TRIBI20MAXN_VpC    3
#define TRIBI20MAXN_CpV    3
#define TRIBI20MAXN_EpC    3
#define TRIBI20MAXN_CpE    2
#define TRIBI20MAXN_EpV    4
#define TRIBI20MAXN_iVpE   1
#define TRIBI20MAXN_nVpoE  3
#define TRIBI20MAXN_nEpoE  2
#define TRIBI20N_E         3

/** refinement descriptor for bisection of a triangle */
class TRefTriBis20Desc : public TRefDesc
{
  public:
    // Constructor
    /** build the descriptor for refining a triangle by bisecting edge 0 */
    TRefTriBis20Desc(TShapeDesc *shape);

    // Methods
};

#endif
