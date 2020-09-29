#ifndef __REFTETRABIS5DESC__
#define __REFTETRABIS5DESC__

#include <RefDesc.h>

#define REFTETRABIS5MAXN_VpC   4
#define REFTETRABIS5MAXN_CpV   2
#define REFTETRABIS5MAXN_EpC   6
#define REFTETRABIS5MAXN_CpE   2
#define REFTETRABIS5MAXN_FpC   4
#define REFTETRABIS5MAXN_CpF   2
#define REFTETRABIS5MAXN_EpV   4
#define REFTETRABIS5MAXN_VpF   3
#define REFTETRABIS5MAXN_FpV   5
#define REFTETRABIS5MAXN_EpF   3
#define REFTETRABIS5MAXN_FpE   3
#define REFTETRABIS5MAXN_nVpoF 4
#define REFTETRABIS5MAXN_oVpoF 3
#define REFTETRABIS5MAXN_iVpE  1
#define REFTETRABIS5MAXN_iEpF  1
#define REFTETRABIS5MAXN_nEpoE 2
#define REFTETRABIS5MAXN_nVpoE 3
#define REFTETRABIS5MAXN_nEpoF 4
#define REFTETRABIS5MAXN_nFpoF 2

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis5Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis5Desc(TShapeDesc *shape);

    // Methods
};

#endif

