#ifndef __REFTETRAQUAD1DESC__
#define __REFTETRAQUAD1DESC__

#include <RefDesc.h>

#define REFTETRAQUAD1MAXN_VpC   4
#define REFTETRAQUAD1MAXN_CpV   4
#define REFTETRAQUAD1MAXN_EpC   6
#define REFTETRAQUAD1MAXN_CpE   3
#define REFTETRAQUAD1MAXN_FpC   4
#define REFTETRAQUAD1MAXN_CpF   2
#define REFTETRAQUAD1MAXN_EpV   6
#define REFTETRAQUAD1MAXN_VpF   3
#define REFTETRAQUAD1MAXN_FpV   9
#define REFTETRAQUAD1MAXN_EpF   3
#define REFTETRAQUAD1MAXN_FpE   4
#define REFTETRAQUAD1MAXN_nVpoF 6
#define REFTETRAQUAD1MAXN_oVpoF 3
#define REFTETRAQUAD1MAXN_iVpE  1
#define REFTETRAQUAD1MAXN_iEpF  3
#define REFTETRAQUAD1MAXN_nEpoE 2
#define REFTETRAQUAD1MAXN_nVpoE 3
#define REFTETRAQUAD1MAXN_nEpoF 6
#define REFTETRAQUAD1MAXN_nFpoF 4

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraQuad1Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraQuad1Desc(TShapeDesc *shape);

    // Methods
};

#endif

