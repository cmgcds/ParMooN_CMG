#ifndef __REFTETRAQUAD0DESC__
#define __REFTETRAQUAD0DESC__

#include <RefDesc.h>

#define REFTETRAQUAD0MAXN_VpC   4
#define REFTETRAQUAD0MAXN_CpV   4
#define REFTETRAQUAD0MAXN_EpC   6
#define REFTETRAQUAD0MAXN_CpE   3
#define REFTETRAQUAD0MAXN_FpC   4
#define REFTETRAQUAD0MAXN_CpF   2
#define REFTETRAQUAD0MAXN_EpV   6
#define REFTETRAQUAD0MAXN_VpF   3
#define REFTETRAQUAD0MAXN_FpV   9
#define REFTETRAQUAD0MAXN_EpF   3
#define REFTETRAQUAD0MAXN_FpE   4
#define REFTETRAQUAD0MAXN_nVpoF 6
#define REFTETRAQUAD0MAXN_oVpoF 3
#define REFTETRAQUAD0MAXN_iVpE  1
#define REFTETRAQUAD0MAXN_iEpF  3
#define REFTETRAQUAD0MAXN_nEpoE 2
#define REFTETRAQUAD0MAXN_nVpoE 3
#define REFTETRAQUAD0MAXN_nEpoF 6
#define REFTETRAQUAD0MAXN_nFpoF 4

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraQuad0Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraQuad0Desc(TShapeDesc *shape);

    // Methods
};

#endif
