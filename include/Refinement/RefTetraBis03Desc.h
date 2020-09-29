
#ifndef __REFTETRABIS03DESC__
#define __REFTETRABIS03DESC__

#include <RefDesc.h>

#define REFTETRABIS03MAXN_VpC   4
#define REFTETRABIS03MAXN_CpV   3
#define REFTETRABIS03MAXN_EpC   6
#define REFTETRABIS03MAXN_CpE   3
#define REFTETRABIS03MAXN_FpC   4
#define REFTETRABIS03MAXN_CpF   2
#define REFTETRABIS03MAXN_EpV   5
#define REFTETRABIS03MAXN_VpF   3
#define REFTETRABIS03MAXN_FpV   7
#define REFTETRABIS03MAXN_EpF   3
#define REFTETRABIS03MAXN_FpE   4
#define REFTETRABIS03MAXN_nVpoF 5
#define REFTETRABIS03MAXN_oVpoF 3
#define REFTETRABIS03MAXN_iVpE  1
#define REFTETRABIS03MAXN_iEpF  2
#define REFTETRABIS03MAXN_nEpoE 2
#define REFTETRABIS03MAXN_nVpoE 3
#define REFTETRABIS03MAXN_nEpoF 5
#define REFTETRABIS03MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis03Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis03Desc(TShapeDesc *shape);

    // Methods
};

#endif
