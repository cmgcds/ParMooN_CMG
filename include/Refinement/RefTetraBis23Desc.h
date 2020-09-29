#ifndef __REFTETRABIS23DESC__
#define __REFTETRABIS23DESC__

#include <RefDesc.h>

#define REFTETRABIS23MAXN_VpC   4
#define REFTETRABIS23MAXN_CpV   3
#define REFTETRABIS23MAXN_EpC   6
#define REFTETRABIS23MAXN_CpE   3
#define REFTETRABIS23MAXN_FpC   4
#define REFTETRABIS23MAXN_CpF   2
#define REFTETRABIS23MAXN_EpV   5
#define REFTETRABIS23MAXN_VpF   3
#define REFTETRABIS23MAXN_FpV   7
#define REFTETRABIS23MAXN_EpF   3
#define REFTETRABIS23MAXN_FpE   4
#define REFTETRABIS23MAXN_nVpoF 5
#define REFTETRABIS23MAXN_oVpoF 3
#define REFTETRABIS23MAXN_iVpE  1
#define REFTETRABIS23MAXN_iEpF  2
#define REFTETRABIS23MAXN_nEpoE 2
#define REFTETRABIS23MAXN_nVpoE 3
#define REFTETRABIS23MAXN_nEpoF 5
#define REFTETRABIS23MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis23Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis23Desc(TShapeDesc *shape);

    // Methods
};

#endif

