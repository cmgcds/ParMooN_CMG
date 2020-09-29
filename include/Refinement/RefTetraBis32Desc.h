#ifndef __REFTETRABIS32DESC__
#define __REFTETRABIS32DESC__

#include <RefDesc.h>

#define REFTETRABIS32MAXN_VpC   4
#define REFTETRABIS32MAXN_CpV   3
#define REFTETRABIS32MAXN_EpC   6
#define REFTETRABIS32MAXN_CpE   3
#define REFTETRABIS32MAXN_FpC   4
#define REFTETRABIS32MAXN_CpF   2
#define REFTETRABIS32MAXN_EpV   5
#define REFTETRABIS32MAXN_VpF   3
#define REFTETRABIS32MAXN_FpV   7
#define REFTETRABIS32MAXN_EpF   3
#define REFTETRABIS32MAXN_FpE   4
#define REFTETRABIS32MAXN_nVpoF 5
#define REFTETRABIS32MAXN_oVpoF 3
#define REFTETRABIS32MAXN_iVpE  1
#define REFTETRABIS32MAXN_iEpF  2
#define REFTETRABIS32MAXN_nEpoE 2
#define REFTETRABIS32MAXN_nVpoE 3
#define REFTETRABIS32MAXN_nEpoF 5
#define REFTETRABIS32MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis32Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis32Desc(TShapeDesc *shape);

    // Methods
};

#endif


