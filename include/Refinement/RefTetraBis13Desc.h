#ifndef __REFTETRABIS13DESC__
#define __REFTETRABIS13DESC__

#include <RefDesc.h>

#define REFTETRABIS13MAXN_VpC   4
#define REFTETRABIS13MAXN_CpV   4
#define REFTETRABIS13MAXN_EpC   6
#define REFTETRABIS13MAXN_CpE   4
#define REFTETRABIS13MAXN_FpC   4
#define REFTETRABIS13MAXN_CpF   2
#define REFTETRABIS13MAXN_EpV   5
#define REFTETRABIS13MAXN_VpF   3
#define REFTETRABIS13MAXN_FpV   9
#define REFTETRABIS13MAXN_EpF   3
#define REFTETRABIS13MAXN_FpE   4
#define REFTETRABIS13MAXN_nVpoF 4
#define REFTETRABIS13MAXN_oVpoF 3
#define REFTETRABIS13MAXN_iVpE  1
#define REFTETRABIS13MAXN_iEpF  1
#define REFTETRABIS13MAXN_nEpoE 2
#define REFTETRABIS13MAXN_nVpoE 3
#define REFTETRABIS13MAXN_nEpoF 4
#define REFTETRABIS13MAXN_nFpoF 2

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis13Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis13Desc(TShapeDesc *shape);

    // Methods
};

#endif

