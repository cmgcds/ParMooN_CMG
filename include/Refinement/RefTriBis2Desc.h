// =======================================================================
// @(#)RefTriBis2Desc.h        1.1 10/30/98
//
// Class:       TRefTriBis2Desc
// Purpose:     refinement descriptor for bisection of a triangle
//              bisection of edge 2
//
// Author:      Volker Behns  18.07.97
//
// =======================================================================

#ifndef __REFTRIBIS2DESC__
#define __REFTRIBIS2DESC__

#include <RefDesc.h>

#define TRIBI2MAXN_VpC    3
#define TRIBI2MAXN_CpV    2
#define TRIBI2MAXN_EpC    3
#define TRIBI2MAXN_CpE    2
#define TRIBI2MAXN_EpV    3
#define TRIBI2MAXN_iVpE   1
#define TRIBI2MAXN_nVpoE  3
#define TRIBI2MAXN_nEpoE  2
#define TRIBI2N_E         3

/** refinement descriptor for bisection of a triangle */
class TRefTriBis2Desc : public TRefDesc
{
  public:
    // Constructor
    /** build the descriptor for refining a triangle by bisecting edge 2 */
    TRefTriBis2Desc(TShapeDesc *shape);

    // Methods
};

#endif
