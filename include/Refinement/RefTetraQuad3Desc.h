#ifndef __REFTETRAQUAD3DESC__
#define __REFTETRAQUAD3DESC__

#include <RefDesc.h>

#define REFTETRAQUAD3MAXN_VpC   4
#define REFTETRAQUAD3MAXN_CpV   4
#define REFTETRAQUAD3MAXN_EpC   6
#define REFTETRAQUAD3MAXN_CpE   3
#define REFTETRAQUAD3MAXN_FpC   4
#define REFTETRAQUAD3MAXN_CpF   2
#define REFTETRAQUAD3MAXN_EpV   6
#define REFTETRAQUAD3MAXN_VpF   3
#define REFTETRAQUAD3MAXN_FpV   9
#define REFTETRAQUAD3MAXN_EpF   3
#define REFTETRAQUAD3MAXN_FpE   4
#define REFTETRAQUAD3MAXN_nVpoF 6
#define REFTETRAQUAD3MAXN_oVpoF 3
#define REFTETRAQUAD3MAXN_iVpE  1
#define REFTETRAQUAD3MAXN_iEpF  3
#define REFTETRAQUAD3MAXN_nEpoE 2
#define REFTETRAQUAD3MAXN_nVpoE 3
#define REFTETRAQUAD3MAXN_nEpoF 6
#define REFTETRAQUAD3MAXN_nFpoF 4

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraQuad3Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraQuad3Desc(TShapeDesc *shape);

    // Methods
};

#endif
