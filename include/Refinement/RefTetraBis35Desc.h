#ifndef __REFTETRABIS35DESC__
#define __REFTETRABIS35DESC__

#include <RefDesc.h>

#define REFTETRABIS35MAXN_VpC   4
#define REFTETRABIS35MAXN_CpV   3
#define REFTETRABIS35MAXN_EpC   6
#define REFTETRABIS35MAXN_CpE   3
#define REFTETRABIS35MAXN_FpC   4
#define REFTETRABIS35MAXN_CpF   2
#define REFTETRABIS35MAXN_EpV   5
#define REFTETRABIS35MAXN_VpF   3
#define REFTETRABIS35MAXN_FpV   7
#define REFTETRABIS35MAXN_EpF   3
#define REFTETRABIS35MAXN_FpE   4
#define REFTETRABIS35MAXN_nVpoF 5
#define REFTETRABIS35MAXN_oVpoF 3
#define REFTETRABIS35MAXN_iVpE  1
#define REFTETRABIS35MAXN_iEpF  2
#define REFTETRABIS35MAXN_nEpoE 2
#define REFTETRABIS35MAXN_nVpoE 3
#define REFTETRABIS35MAXN_nEpoF 5
#define REFTETRABIS35MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis35Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis35Desc(TShapeDesc *shape);

    // Methods
};

#endif
