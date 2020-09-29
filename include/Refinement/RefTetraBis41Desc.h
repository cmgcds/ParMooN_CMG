#ifndef __REFTETRABIS41DESC__
#define __REFTETRABIS41DESC__

#include <RefDesc.h>

#define REFTETRABIS41MAXN_VpC   4
#define REFTETRABIS41MAXN_CpV   3
#define REFTETRABIS41MAXN_EpC   6
#define REFTETRABIS41MAXN_CpE   3
#define REFTETRABIS41MAXN_FpC   4
#define REFTETRABIS41MAXN_CpF   2
#define REFTETRABIS41MAXN_EpV   5
#define REFTETRABIS41MAXN_VpF   3
#define REFTETRABIS41MAXN_FpV   7
#define REFTETRABIS41MAXN_EpF   3
#define REFTETRABIS41MAXN_FpE   4
#define REFTETRABIS41MAXN_nVpoF 5
#define REFTETRABIS41MAXN_oVpoF 3
#define REFTETRABIS41MAXN_iVpE  1
#define REFTETRABIS41MAXN_iEpF  2
#define REFTETRABIS41MAXN_nEpoE 2
#define REFTETRABIS41MAXN_nVpoE 3
#define REFTETRABIS41MAXN_nEpoF 5
#define REFTETRABIS41MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis41Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis41Desc(TShapeDesc *shape);

    // Methods
};

#endif

