#ifndef __REFTETRABIS2DESC__
#define __REFTETRABIS2DESC__

#include <RefDesc.h>

#define REFTETRABIS2MAXN_VpC   4
#define REFTETRABIS2MAXN_CpV   2
#define REFTETRABIS2MAXN_EpC   6
#define REFTETRABIS2MAXN_CpE   2
#define REFTETRABIS2MAXN_FpC   4
#define REFTETRABIS2MAXN_CpF   2
#define REFTETRABIS2MAXN_EpV   4
#define REFTETRABIS2MAXN_VpF   3
#define REFTETRABIS2MAXN_FpV   5
#define REFTETRABIS2MAXN_EpF   3
#define REFTETRABIS2MAXN_FpE   3
#define REFTETRABIS2MAXN_nVpoF 4
#define REFTETRABIS2MAXN_oVpoF 3
#define REFTETRABIS2MAXN_iVpE  1
#define REFTETRABIS2MAXN_iEpF  1
#define REFTETRABIS2MAXN_nEpoE 2
#define REFTETRABIS2MAXN_nVpoE 3
#define REFTETRABIS2MAXN_nEpoF 4
#define REFTETRABIS2MAXN_nFpoF 2

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis2Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis2Desc(TShapeDesc *shape);

    // Methods
};

#endif

