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

#ifndef __REFTRIBIS10DESC__
#define __REFTRIBIS10DESC__

#include <RefDesc.h>

#define TRIBI10MAXN_VpC    3
#define TRIBI10MAXN_CpV    3
#define TRIBI10MAXN_EpC    3
#define TRIBI10MAXN_CpE    2
#define TRIBI10MAXN_EpV    4
#define TRIBI10MAXN_iVpE   1
#define TRIBI10MAXN_nVpoE  3
#define TRIBI10MAXN_nEpoE  2
#define TRIBI10N_E         3

/** refinement descriptor for bisection of a triangle */
class TRefTriBis10Desc : public TRefDesc
{
  public:
    // Constructor
    /** build the descriptor for refining a triangle by bisecting edge 0 */
    TRefTriBis10Desc(TShapeDesc *shape);

    // Methods
};

#endif
