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

#ifndef __REFTRIBIS02DESC__
#define __REFTRIBIS02DESC__

#include <RefDesc.h>

#define TRIBI02MAXN_VpC    3
#define TRIBI02MAXN_CpV    3
#define TRIBI02MAXN_EpC    3
#define TRIBI02MAXN_CpE    2
#define TRIBI02MAXN_EpV    4
#define TRIBI02MAXN_iVpE   1
#define TRIBI02MAXN_nVpoE  3
#define TRIBI02MAXN_nEpoE  2
#define TRIBI02N_E         3

/** refinement descriptor for bisection of a triangle */
class TRefTriBis02Desc : public TRefDesc
{
  public:
    // Constructor
    /** build the descriptor for refining a triangle by bisecting edge 0 */
    TRefTriBis02Desc(TShapeDesc *shape);

    // Methods
};

#endif
