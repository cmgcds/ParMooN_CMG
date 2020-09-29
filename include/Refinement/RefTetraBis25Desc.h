#ifndef __REFTETRABIS25DESC__
#define __REFTETRABIS25DESC__

#include <RefDesc.h>

#define REFTETRABIS25MAXN_VpC   4
#define REFTETRABIS25MAXN_CpV   3
#define REFTETRABIS25MAXN_EpC   6
#define REFTETRABIS25MAXN_CpE   3
#define REFTETRABIS25MAXN_FpC   4
#define REFTETRABIS25MAXN_CpF   2
#define REFTETRABIS25MAXN_EpV   5
#define REFTETRABIS25MAXN_VpF   3
#define REFTETRABIS25MAXN_FpV   7
#define REFTETRABIS25MAXN_EpF   3
#define REFTETRABIS25MAXN_FpE   4
#define REFTETRABIS25MAXN_nVpoF 5
#define REFTETRABIS25MAXN_oVpoF 3
#define REFTETRABIS25MAXN_iVpE  1
#define REFTETRABIS25MAXN_iEpF  2
#define REFTETRABIS25MAXN_nEpoE 2
#define REFTETRABIS25MAXN_nVpoE 3
#define REFTETRABIS25MAXN_nEpoF 5
#define REFTETRABIS25MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis25Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis25Desc(TShapeDesc *shape);

    // Methods
};

#endif

