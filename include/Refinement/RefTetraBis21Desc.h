#ifndef __REFTETRABIS21DESC__
#define __REFTETRABIS21DESC__

#include <RefDesc.h>

#define REFTETRABIS21MAXN_VpC   4
#define REFTETRABIS21MAXN_CpV   3
#define REFTETRABIS21MAXN_EpC   6
#define REFTETRABIS21MAXN_CpE   3
#define REFTETRABIS21MAXN_FpC   4
#define REFTETRABIS21MAXN_CpF   2
#define REFTETRABIS21MAXN_EpV   5
#define REFTETRABIS21MAXN_VpF   3
#define REFTETRABIS21MAXN_FpV   7
#define REFTETRABIS21MAXN_EpF   3
#define REFTETRABIS21MAXN_FpE   4
#define REFTETRABIS21MAXN_nVpoF 5
#define REFTETRABIS21MAXN_oVpoF 3
#define REFTETRABIS21MAXN_iVpE  1
#define REFTETRABIS21MAXN_iEpF  2
#define REFTETRABIS21MAXN_nEpoE 2
#define REFTETRABIS21MAXN_nVpoE 3
#define REFTETRABIS21MAXN_nEpoF 5
#define REFTETRABIS21MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis21Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis21Desc(TShapeDesc *shape);

    // Methods
};

#endif

