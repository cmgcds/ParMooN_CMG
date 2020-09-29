#ifndef __REFTETRABIS4DESC__
#define __REFTETRABIS4DESC__

#include <RefDesc.h>

#define REFTETRABIS4MAXN_VpC   4
#define REFTETRABIS4MAXN_CpV   2
#define REFTETRABIS4MAXN_EpC   6
#define REFTETRABIS4MAXN_CpE   2
#define REFTETRABIS4MAXN_FpC   4
#define REFTETRABIS4MAXN_CpF   2
#define REFTETRABIS4MAXN_EpV   4
#define REFTETRABIS4MAXN_VpF   3
#define REFTETRABIS4MAXN_FpV   5
#define REFTETRABIS4MAXN_EpF   3
#define REFTETRABIS4MAXN_FpE   3
#define REFTETRABIS4MAXN_nVpoF 4
#define REFTETRABIS4MAXN_oVpoF 3
#define REFTETRABIS4MAXN_iVpE  1
#define REFTETRABIS4MAXN_iEpF  1
#define REFTETRABIS4MAXN_nEpoE 2
#define REFTETRABIS4MAXN_nVpoE 3
#define REFTETRABIS4MAXN_nEpoF 4
#define REFTETRABIS4MAXN_nFpoF 2

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis4Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis4Desc(TShapeDesc *shape);

    // Methods
};

#endif

