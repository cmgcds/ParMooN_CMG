#ifndef __REFTETRABIS02DESC__
#define __REFTETRABIS02DESC__

#include <RefDesc.h>

#define REFTETRABIS02MAXN_VpC   4
#define REFTETRABIS02MAXN_CpV   3
#define REFTETRABIS02MAXN_EpC   6
#define REFTETRABIS02MAXN_CpE   3
#define REFTETRABIS02MAXN_FpC   4
#define REFTETRABIS02MAXN_CpF   2
#define REFTETRABIS02MAXN_EpV   5
#define REFTETRABIS02MAXN_VpF   3
#define REFTETRABIS02MAXN_FpV   7
#define REFTETRABIS02MAXN_EpF   3
#define REFTETRABIS02MAXN_FpE   4
#define REFTETRABIS02MAXN_nVpoF 5
#define REFTETRABIS02MAXN_oVpoF 3
#define REFTETRABIS02MAXN_iVpE  1
#define REFTETRABIS02MAXN_iEpF  2
#define REFTETRABIS02MAXN_nEpoE 2
#define REFTETRABIS02MAXN_nVpoE 3
#define REFTETRABIS02MAXN_nEpoF 5
#define REFTETRABIS02MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis02Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis02Desc(TShapeDesc *shape);

    // Methods
};

#endif

