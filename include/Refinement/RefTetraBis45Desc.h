#ifndef __REFTETRABIS45DESC__
#define __REFTETRABIS45DESC__

#include <RefDesc.h>

#define REFTETRABIS45MAXN_VpC   4
#define REFTETRABIS45MAXN_CpV   3
#define REFTETRABIS45MAXN_EpC   6
#define REFTETRABIS45MAXN_CpE   3
#define REFTETRABIS45MAXN_FpC   4
#define REFTETRABIS45MAXN_CpF   2
#define REFTETRABIS45MAXN_EpV   5
#define REFTETRABIS45MAXN_VpF   3
#define REFTETRABIS45MAXN_FpV   7
#define REFTETRABIS45MAXN_EpF   3
#define REFTETRABIS45MAXN_FpE   4
#define REFTETRABIS45MAXN_nVpoF 5
#define REFTETRABIS45MAXN_oVpoF 3
#define REFTETRABIS45MAXN_iVpE  1
#define REFTETRABIS45MAXN_iEpF  2
#define REFTETRABIS45MAXN_nEpoE 2
#define REFTETRABIS45MAXN_nVpoE 3
#define REFTETRABIS45MAXN_nEpoF 5
#define REFTETRABIS45MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis45Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis45Desc(TShapeDesc *shape);

    // Methods
};

#endif

