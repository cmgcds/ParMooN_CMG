#ifndef __REFTETRABIS20DESC__
#define __REFTETRABIS20DESC__

#include <RefDesc.h>

#define REFTETRABIS20MAXN_VpC   4
#define REFTETRABIS20MAXN_CpV   3
#define REFTETRABIS20MAXN_EpC   6
#define REFTETRABIS20MAXN_CpE   3
#define REFTETRABIS20MAXN_FpC   4
#define REFTETRABIS20MAXN_CpF   2
#define REFTETRABIS20MAXN_EpV   5
#define REFTETRABIS20MAXN_VpF   3
#define REFTETRABIS20MAXN_FpV   7
#define REFTETRABIS20MAXN_EpF   3
#define REFTETRABIS20MAXN_FpE   4
#define REFTETRABIS20MAXN_nVpoF 5
#define REFTETRABIS20MAXN_oVpoF 3
#define REFTETRABIS20MAXN_iVpE  1
#define REFTETRABIS20MAXN_iEpF  2
#define REFTETRABIS20MAXN_nEpoE 2
#define REFTETRABIS20MAXN_nVpoE 3
#define REFTETRABIS20MAXN_nEpoF 5
#define REFTETRABIS20MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis20Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis20Desc(TShapeDesc *shape);

    // Methods
};

#endif
