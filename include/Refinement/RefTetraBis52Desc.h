#ifndef __REFTETRABIS52DESC__
#define __REFTETRABIS52DESC__

#include <RefDesc.h>

#define REFTETRABIS52MAXN_VpC   4
#define REFTETRABIS52MAXN_CpV   3
#define REFTETRABIS52MAXN_EpC   6
#define REFTETRABIS52MAXN_CpE   3
#define REFTETRABIS52MAXN_FpC   4
#define REFTETRABIS52MAXN_CpF   2
#define REFTETRABIS52MAXN_EpV   5
#define REFTETRABIS52MAXN_VpF   3
#define REFTETRABIS52MAXN_FpV   7
#define REFTETRABIS52MAXN_EpF   3
#define REFTETRABIS52MAXN_FpE   4
#define REFTETRABIS52MAXN_nVpoF 5
#define REFTETRABIS52MAXN_oVpoF 3
#define REFTETRABIS52MAXN_iVpE  1
#define REFTETRABIS52MAXN_iEpF  2
#define REFTETRABIS52MAXN_nEpoE 2
#define REFTETRABIS52MAXN_nVpoE 3
#define REFTETRABIS52MAXN_nEpoF 5
#define REFTETRABIS52MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis52Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis52Desc(TShapeDesc *shape);

    // Methods
};

#endif

