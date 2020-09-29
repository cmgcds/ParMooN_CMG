#ifndef __REFTETRAQUAD2DESC__
#define __REFTETRAQUAD2DESC__

#include <RefDesc.h>

#define REFTETRAQUAD2MAXN_VpC   4
#define REFTETRAQUAD2MAXN_CpV   4
#define REFTETRAQUAD2MAXN_EpC   6
#define REFTETRAQUAD2MAXN_CpE   3
#define REFTETRAQUAD2MAXN_FpC   4
#define REFTETRAQUAD2MAXN_CpF   2
#define REFTETRAQUAD2MAXN_EpV   6
#define REFTETRAQUAD2MAXN_VpF   3
#define REFTETRAQUAD2MAXN_FpV   9
#define REFTETRAQUAD2MAXN_EpF   3
#define REFTETRAQUAD2MAXN_FpE   4
#define REFTETRAQUAD2MAXN_nVpoF 6
#define REFTETRAQUAD2MAXN_oVpoF 3
#define REFTETRAQUAD2MAXN_iVpE  1
#define REFTETRAQUAD2MAXN_iEpF  3
#define REFTETRAQUAD2MAXN_nEpoE 2
#define REFTETRAQUAD2MAXN_nVpoE 3
#define REFTETRAQUAD2MAXN_nEpoF 6
#define REFTETRAQUAD2MAXN_nFpoF 4

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraQuad2Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraQuad2Desc(TShapeDesc *shape);

    // Methods
};

#endif
