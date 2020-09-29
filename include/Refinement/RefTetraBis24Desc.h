#ifndef __REFTETRABIS24DESC__
#define __REFTETRABIS24DESC__

#include <RefDesc.h>

#define REFTETRABIS24MAXN_VpC   4
#define REFTETRABIS24MAXN_CpV   4
#define REFTETRABIS24MAXN_EpC   6
#define REFTETRABIS24MAXN_CpE   4
#define REFTETRABIS24MAXN_FpC   4
#define REFTETRABIS24MAXN_CpF   2
#define REFTETRABIS24MAXN_EpV   5
#define REFTETRABIS24MAXN_VpF   3
#define REFTETRABIS24MAXN_FpV   9
#define REFTETRABIS24MAXN_EpF   3
#define REFTETRABIS24MAXN_FpE   4
#define REFTETRABIS24MAXN_nVpoF 4
#define REFTETRABIS24MAXN_oVpoF 3
#define REFTETRABIS24MAXN_iVpE  1
#define REFTETRABIS24MAXN_iEpF  1
#define REFTETRABIS24MAXN_nEpoE 2
#define REFTETRABIS24MAXN_nVpoE 3
#define REFTETRABIS24MAXN_nEpoF 4
#define REFTETRABIS24MAXN_nFpoF 2

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis24Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis24Desc(TShapeDesc *shape);

    // Methods
};

#endif


