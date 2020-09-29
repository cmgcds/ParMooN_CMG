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

#ifndef __REFTRIBIS21DESC__
#define __REFTRIBIS21DESC__

#include <RefDesc.h>

#define TRIBI21MAXN_VpC    3
#define TRIBI21MAXN_CpV    3
#define TRIBI21MAXN_EpC    3
#define TRIBI21MAXN_CpE    2
#define TRIBI21MAXN_EpV    4
#define TRIBI21MAXN_iVpE   1
#define TRIBI21MAXN_nVpoE  3
#define TRIBI21MAXN_nEpoE  2
#define TRIBI21N_E         3

/** refinement descriptor for bisection of a triangle */
class TRefTriBis21Desc : public TRefDesc
{
  public:
    // Constructor
    /** build the descriptor for refining a triangle by bisecting edge 0 */
    TRefTriBis21Desc(TShapeDesc *shape);

    // Methods
};

#endif
