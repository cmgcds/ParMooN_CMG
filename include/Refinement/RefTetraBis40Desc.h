#ifndef __REFTETRABIS40DESC__
#define __REFTETRABIS40DESC__

#include <RefDesc.h>

#define REFTETRABIS40MAXN_VpC   4
#define REFTETRABIS40MAXN_CpV   3
#define REFTETRABIS40MAXN_EpC   6
#define REFTETRABIS40MAXN_CpE   3
#define REFTETRABIS40MAXN_FpC   4
#define REFTETRABIS40MAXN_CpF   2
#define REFTETRABIS40MAXN_EpV   5
#define REFTETRABIS40MAXN_VpF   3
#define REFTETRABIS40MAXN_FpV   7
#define REFTETRABIS40MAXN_EpF   3
#define REFTETRABIS40MAXN_FpE   4
#define REFTETRABIS40MAXN_nVpoF 5
#define REFTETRABIS40MAXN_oVpoF 3
#define REFTETRABIS40MAXN_iVpE  1
#define REFTETRABIS40MAXN_iEpF  2
#define REFTETRABIS40MAXN_nEpoE 2
#define REFTETRABIS40MAXN_nVpoE 3
#define REFTETRABIS40MAXN_nEpoF 5
#define REFTETRABIS40MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis40Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis40Desc(TShapeDesc *shape);

    // Methods
};

#endif
