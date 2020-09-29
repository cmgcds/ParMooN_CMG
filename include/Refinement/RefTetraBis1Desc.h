#ifndef __REFTETRABIS1DESC__
#define __REFTETRABIS1DESC__

#include <RefDesc.h>

#define REFTETRABIS1MAXN_VpC   4
#define REFTETRABIS1MAXN_CpV   2
#define REFTETRABIS1MAXN_EpC   6
#define REFTETRABIS1MAXN_CpE   2
#define REFTETRABIS1MAXN_FpC   4
#define REFTETRABIS1MAXN_CpF   2
#define REFTETRABIS1MAXN_EpV   4
#define REFTETRABIS1MAXN_VpF   3
#define REFTETRABIS1MAXN_FpV   5
#define REFTETRABIS1MAXN_EpF   3
#define REFTETRABIS1MAXN_FpE   3
#define REFTETRABIS1MAXN_nVpoF 4
#define REFTETRABIS1MAXN_oVpoF 3
#define REFTETRABIS1MAXN_iVpE  1
#define REFTETRABIS1MAXN_iEpF  1
#define REFTETRABIS1MAXN_nEpoE 2
#define REFTETRABIS1MAXN_nVpoE 3
#define REFTETRABIS1MAXN_nEpoF 5
#define REFTETRABIS1MAXN_nFpoF 2

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis1Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis1Desc(TShapeDesc *shape);

    // Methods
};

#endif
