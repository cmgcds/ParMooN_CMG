#ifndef __REFTETRABIS12DESC__
#define __REFTETRABIS12DESC__

#include <RefDesc.h>

#define REFTETRABIS12MAXN_VpC   4
#define REFTETRABIS12MAXN_CpV   3
#define REFTETRABIS12MAXN_EpC   6
#define REFTETRABIS12MAXN_CpE   3
#define REFTETRABIS12MAXN_FpC   4
#define REFTETRABIS12MAXN_CpF   2
#define REFTETRABIS12MAXN_EpV   5
#define REFTETRABIS12MAXN_VpF   3
#define REFTETRABIS12MAXN_FpV   7
#define REFTETRABIS12MAXN_EpF   3
#define REFTETRABIS12MAXN_FpE   4
#define REFTETRABIS12MAXN_nVpoF 5
#define REFTETRABIS12MAXN_oVpoF 3
#define REFTETRABIS12MAXN_iVpE  1
#define REFTETRABIS12MAXN_iEpF  2
#define REFTETRABIS12MAXN_nEpoE 2
#define REFTETRABIS12MAXN_nVpoE 3
#define REFTETRABIS12MAXN_nEpoF 5
#define REFTETRABIS12MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis12Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis12Desc(TShapeDesc *shape);

    // Methods
};

#endif

