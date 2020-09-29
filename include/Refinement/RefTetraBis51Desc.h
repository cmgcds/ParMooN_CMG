#ifndef __REFTETRABIS51DESC__
#define __REFTETRABIS51DESC__

#include <RefDesc.h>

#define REFTETRABIS51MAXN_VpC   4
#define REFTETRABIS51MAXN_CpV   3
#define REFTETRABIS51MAXN_EpC   6
#define REFTETRABIS51MAXN_CpE   3
#define REFTETRABIS51MAXN_FpC   4
#define REFTETRABIS51MAXN_CpF   2
#define REFTETRABIS51MAXN_EpV   5
#define REFTETRABIS51MAXN_VpF   3
#define REFTETRABIS51MAXN_FpV   7
#define REFTETRABIS51MAXN_EpF   3
#define REFTETRABIS51MAXN_FpE   4
#define REFTETRABIS51MAXN_nVpoF 5
#define REFTETRABIS51MAXN_oVpoF 3
#define REFTETRABIS51MAXN_iVpE  1
#define REFTETRABIS51MAXN_iEpF  2
#define REFTETRABIS51MAXN_nEpoE 2
#define REFTETRABIS51MAXN_nVpoE 3
#define REFTETRABIS51MAXN_nEpoF 5
#define REFTETRABIS51MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis51Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis51Desc(TShapeDesc *shape);

    // Methods
};

#endif

