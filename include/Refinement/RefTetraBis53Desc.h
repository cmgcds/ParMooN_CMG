#ifndef __REFTETRABIS53DESC__
#define __REFTETRABIS53DESC__

#include <RefDesc.h>

#define REFTETRABIS53MAXN_VpC   4
#define REFTETRABIS53MAXN_CpV   3
#define REFTETRABIS53MAXN_EpC   6
#define REFTETRABIS53MAXN_CpE   3
#define REFTETRABIS53MAXN_FpC   4
#define REFTETRABIS53MAXN_CpF   2
#define REFTETRABIS53MAXN_EpV   5
#define REFTETRABIS53MAXN_VpF   3
#define REFTETRABIS53MAXN_FpV   7
#define REFTETRABIS53MAXN_EpF   3
#define REFTETRABIS53MAXN_FpE   4
#define REFTETRABIS53MAXN_nVpoF 5
#define REFTETRABIS53MAXN_oVpoF 3
#define REFTETRABIS53MAXN_iVpE  1
#define REFTETRABIS53MAXN_iEpF  2
#define REFTETRABIS53MAXN_nEpoE 2
#define REFTETRABIS53MAXN_nVpoE 3
#define REFTETRABIS53MAXN_nEpoF 5
#define REFTETRABIS53MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis53Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis53Desc(TShapeDesc *shape);

    // Methods
};

#endif

