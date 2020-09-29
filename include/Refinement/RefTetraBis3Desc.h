#ifndef __REFTETRABIS3DESC__
#define __REFTETRABIS3DESC__

#include <RefDesc.h>

#define REFTETRABIS3MAXN_VpC   4
#define REFTETRABIS3MAXN_CpV   2
#define REFTETRABIS3MAXN_EpC   6
#define REFTETRABIS3MAXN_CpE   2
#define REFTETRABIS3MAXN_FpC   4
#define REFTETRABIS3MAXN_CpF   2
#define REFTETRABIS3MAXN_EpV   4
#define REFTETRABIS3MAXN_VpF   3
#define REFTETRABIS3MAXN_FpV   5
#define REFTETRABIS3MAXN_EpF   3
#define REFTETRABIS3MAXN_FpE   3
#define REFTETRABIS3MAXN_nVpoF 4
#define REFTETRABIS3MAXN_oVpoF 3
#define REFTETRABIS3MAXN_iVpE  1
#define REFTETRABIS3MAXN_iEpF  1
#define REFTETRABIS3MAXN_nEpoE 2
#define REFTETRABIS3MAXN_nVpoE 3
#define REFTETRABIS3MAXN_nEpoF 4
#define REFTETRABIS3MAXN_nFpoF 2

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis3Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis3Desc(TShapeDesc *shape);

    // Methods
};

#endif

