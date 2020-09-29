#ifndef __REFTETRABIS05DESC__
#define __REFTETRABIS05DESC__

#include <RefDesc.h>

#define REFTETRABIS05MAXN_VpC   4
#define REFTETRABIS05MAXN_CpV   4
#define REFTETRABIS05MAXN_EpC   6
#define REFTETRABIS05MAXN_CpE   4
#define REFTETRABIS05MAXN_FpC   4
#define REFTETRABIS05MAXN_CpF   2
#define REFTETRABIS05MAXN_EpV   5
#define REFTETRABIS05MAXN_VpF   3
#define REFTETRABIS05MAXN_FpV   9
#define REFTETRABIS05MAXN_EpF   3
#define REFTETRABIS05MAXN_FpE   4
#define REFTETRABIS05MAXN_nVpoF 4
#define REFTETRABIS05MAXN_oVpoF 3
#define REFTETRABIS05MAXN_iVpE  1
#define REFTETRABIS05MAXN_iEpF  1
#define REFTETRABIS05MAXN_nEpoE 2
#define REFTETRABIS05MAXN_nVpoE 3
#define REFTETRABIS05MAXN_nEpoF 4
#define REFTETRABIS05MAXN_nFpoF 2

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis05Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis05Desc(TShapeDesc *shape);

    // Methods
};

#endif
