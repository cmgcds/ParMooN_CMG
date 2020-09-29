#ifndef __REFTETRABIS30DESC__
#define __REFTETRABIS30DESC__

#include <RefDesc.h>

#define REFTETRABIS30MAXN_VpC   4
#define REFTETRABIS30MAXN_CpV   3
#define REFTETRABIS30MAXN_EpC   6
#define REFTETRABIS30MAXN_CpE   3
#define REFTETRABIS30MAXN_FpC   4
#define REFTETRABIS30MAXN_CpF   2
#define REFTETRABIS30MAXN_EpV   5
#define REFTETRABIS30MAXN_VpF   3
#define REFTETRABIS30MAXN_FpV   7
#define REFTETRABIS30MAXN_EpF   3
#define REFTETRABIS30MAXN_FpE   4
#define REFTETRABIS30MAXN_nVpoF 5
#define REFTETRABIS30MAXN_oVpoF 3
#define REFTETRABIS30MAXN_iVpE  1
#define REFTETRABIS30MAXN_iEpF  2
#define REFTETRABIS30MAXN_nEpoE 2
#define REFTETRABIS30MAXN_nVpoE 3
#define REFTETRABIS30MAXN_nEpoF 5
#define REFTETRABIS30MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis30Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis30Desc(TShapeDesc *shape);

    // Methods
};

#endif


