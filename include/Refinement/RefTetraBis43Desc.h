#ifndef __REFTETRABIS43DESC__
#define __REFTETRABIS43DESC__

#include <RefDesc.h>

#define REFTETRABIS43MAXN_VpC   4
#define REFTETRABIS43MAXN_CpV   3
#define REFTETRABIS43MAXN_EpC   6
#define REFTETRABIS43MAXN_CpE   3
#define REFTETRABIS43MAXN_FpC   4
#define REFTETRABIS43MAXN_CpF   2
#define REFTETRABIS43MAXN_EpV   5
#define REFTETRABIS43MAXN_VpF   3
#define REFTETRABIS43MAXN_FpV   7
#define REFTETRABIS43MAXN_EpF   3
#define REFTETRABIS43MAXN_FpE   4
#define REFTETRABIS43MAXN_nVpoF 5
#define REFTETRABIS43MAXN_oVpoF 3
#define REFTETRABIS43MAXN_iVpE  1
#define REFTETRABIS43MAXN_iEpF  2
#define REFTETRABIS43MAXN_nEpoE 2
#define REFTETRABIS43MAXN_nVpoE 3
#define REFTETRABIS43MAXN_nEpoF 5
#define REFTETRABIS43MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis43Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis43Desc(TShapeDesc *shape);

    // Methods
};

#endif
