#ifndef __REFTETRABIS01DESC__
#define __REFTETRABIS01DESC__

#include <RefDesc.h>

#define REFTETRABIS01MAXN_VpC   4
#define REFTETRABIS01MAXN_CpV   3
#define REFTETRABIS01MAXN_EpC   6
#define REFTETRABIS01MAXN_CpE   3
#define REFTETRABIS01MAXN_FpC   4
#define REFTETRABIS01MAXN_CpF   2
#define REFTETRABIS01MAXN_EpV   5
#define REFTETRABIS01MAXN_VpF   3
#define REFTETRABIS01MAXN_FpV   7
#define REFTETRABIS01MAXN_EpF   3
#define REFTETRABIS01MAXN_FpE   4
#define REFTETRABIS01MAXN_nVpoF 5
#define REFTETRABIS01MAXN_oVpoF 3
#define REFTETRABIS01MAXN_iVpE  1
#define REFTETRABIS01MAXN_iEpF  2
#define REFTETRABIS01MAXN_nEpoE 2
#define REFTETRABIS01MAXN_nVpoE 3
#define REFTETRABIS01MAXN_nEpoF 5
#define REFTETRABIS01MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis01Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis01Desc(TShapeDesc *shape);

    // Methods
};

#endif
