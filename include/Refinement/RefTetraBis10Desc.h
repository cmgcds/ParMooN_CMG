#ifndef __REFTETRABIS10DESC__
#define __REFTETRABIS10DESC__

#include <RefDesc.h>

#define REFTETRABIS10MAXN_VpC   4
#define REFTETRABIS10MAXN_CpV   3
#define REFTETRABIS10MAXN_EpC   6
#define REFTETRABIS10MAXN_CpE   3
#define REFTETRABIS10MAXN_FpC   4
#define REFTETRABIS10MAXN_CpF   2
#define REFTETRABIS10MAXN_EpV   5
#define REFTETRABIS10MAXN_VpF   3
#define REFTETRABIS10MAXN_FpV   7
#define REFTETRABIS10MAXN_EpF   3
#define REFTETRABIS10MAXN_FpE   4
#define REFTETRABIS10MAXN_nVpoF 5
#define REFTETRABIS10MAXN_oVpoF 3
#define REFTETRABIS10MAXN_iVpE  1
#define REFTETRABIS10MAXN_iEpF  2
#define REFTETRABIS10MAXN_nEpoE 2
#define REFTETRABIS10MAXN_nVpoE 3
#define REFTETRABIS10MAXN_nEpoF 5
#define REFTETRABIS10MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis10Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis10Desc(TShapeDesc *shape);

    // Methods
};

#endif
