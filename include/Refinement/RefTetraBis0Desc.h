#ifndef __REFTETRABIS0DESC__
#define __REFTETRABIS0DESC__

#include <RefDesc.h>

#define REFTETRABIS0MAXN_VpC   4
#define REFTETRABIS0MAXN_CpV   2
#define REFTETRABIS0MAXN_EpC   6
#define REFTETRABIS0MAXN_CpE   2
#define REFTETRABIS0MAXN_FpC   4
#define REFTETRABIS0MAXN_CpF   2
#define REFTETRABIS0MAXN_EpV   4
#define REFTETRABIS0MAXN_VpF   3
#define REFTETRABIS0MAXN_FpV   5
#define REFTETRABIS0MAXN_EpF   3
#define REFTETRABIS0MAXN_FpE   3
#define REFTETRABIS0MAXN_nVpoF 4
#define REFTETRABIS0MAXN_oVpoF 3
#define REFTETRABIS0MAXN_iVpE  1
#define REFTETRABIS0MAXN_iEpF  1
#define REFTETRABIS0MAXN_nEpoE 2
#define REFTETRABIS0MAXN_nVpoE 3
#define REFTETRABIS0MAXN_nEpoF 4
#define REFTETRABIS0MAXN_nFpoF 2

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis0Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis0Desc(TShapeDesc *shape);

    // Methods
};

#endif
