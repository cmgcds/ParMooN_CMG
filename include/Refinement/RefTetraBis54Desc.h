#ifndef __REFTETRABIS54DESC__
#define __REFTETRABIS54DESC__

#include <RefDesc.h>

#define REFTETRABIS54MAXN_VpC   4
#define REFTETRABIS54MAXN_CpV   3
#define REFTETRABIS54MAXN_EpC   6
#define REFTETRABIS54MAXN_CpE   3
#define REFTETRABIS54MAXN_FpC   4
#define REFTETRABIS54MAXN_CpF   2
#define REFTETRABIS54MAXN_EpV   5
#define REFTETRABIS54MAXN_VpF   3
#define REFTETRABIS54MAXN_FpV   7
#define REFTETRABIS54MAXN_EpF   3
#define REFTETRABIS54MAXN_FpE   4
#define REFTETRABIS54MAXN_nVpoF 5
#define REFTETRABIS54MAXN_oVpoF 3
#define REFTETRABIS54MAXN_iVpE  1
#define REFTETRABIS54MAXN_iEpF  2
#define REFTETRABIS54MAXN_nEpoE 2
#define REFTETRABIS54MAXN_nVpoE 3
#define REFTETRABIS54MAXN_nEpoF 5
#define REFTETRABIS54MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis54Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis54Desc(TShapeDesc *shape);

    // Methods
};

#endif

