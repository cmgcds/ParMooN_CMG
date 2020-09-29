#ifndef __REFTETRABIS34DESC__
#define __REFTETRABIS34DESC__

#include <RefDesc.h>

#define REFTETRABIS34MAXN_VpC   4
#define REFTETRABIS34MAXN_CpV   3
#define REFTETRABIS34MAXN_EpC   6
#define REFTETRABIS34MAXN_CpE   3
#define REFTETRABIS34MAXN_FpC   4
#define REFTETRABIS34MAXN_CpF   2
#define REFTETRABIS34MAXN_EpV   5
#define REFTETRABIS34MAXN_VpF   3
#define REFTETRABIS34MAXN_FpV   7
#define REFTETRABIS34MAXN_EpF   3
#define REFTETRABIS34MAXN_FpE   4
#define REFTETRABIS34MAXN_nVpoF 5
#define REFTETRABIS34MAXN_oVpoF 3
#define REFTETRABIS34MAXN_iVpE  1
#define REFTETRABIS34MAXN_iEpF  2
#define REFTETRABIS34MAXN_nEpoE 2
#define REFTETRABIS34MAXN_nVpoE 3
#define REFTETRABIS34MAXN_nEpoF 5
#define REFTETRABIS34MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis34Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis34Desc(TShapeDesc *shape);

    // Methods
};

#endif


