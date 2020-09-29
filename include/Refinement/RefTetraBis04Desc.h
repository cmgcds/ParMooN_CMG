#ifndef __REFTETRABIS04DESC__
#define __REFTETRABIS04DESC__

#include <RefDesc.h>

#define REFTETRABIS04MAXN_VpC   4
#define REFTETRABIS04MAXN_CpV   3
#define REFTETRABIS04MAXN_EpC   6
#define REFTETRABIS04MAXN_CpE   3
#define REFTETRABIS04MAXN_FpC   4
#define REFTETRABIS04MAXN_CpF   2
#define REFTETRABIS04MAXN_EpV   5
#define REFTETRABIS04MAXN_VpF   3
#define REFTETRABIS04MAXN_FpV   7
#define REFTETRABIS04MAXN_EpF   3
#define REFTETRABIS04MAXN_FpE   4
#define REFTETRABIS04MAXN_nVpoF 5
#define REFTETRABIS04MAXN_oVpoF 3
#define REFTETRABIS04MAXN_iVpE  1
#define REFTETRABIS04MAXN_iEpF  2
#define REFTETRABIS04MAXN_nEpoE 2
#define REFTETRABIS04MAXN_nVpoE 3
#define REFTETRABIS04MAXN_nEpoF 5
#define REFTETRABIS04MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis04Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis04Desc(TShapeDesc *shape);

    // Methods
};

#endif

