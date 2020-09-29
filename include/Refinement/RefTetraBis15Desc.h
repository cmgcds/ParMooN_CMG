#ifndef __REFTETRABIS15DESC__
#define __REFTETRABIS15DESC__

#include <RefDesc.h>

#define REFTETRABIS15MAXN_VpC   4
#define REFTETRABIS15MAXN_CpV   3
#define REFTETRABIS15MAXN_EpC   6
#define REFTETRABIS15MAXN_CpE   3
#define REFTETRABIS15MAXN_FpC   4
#define REFTETRABIS15MAXN_CpF   2
#define REFTETRABIS15MAXN_EpV   5
#define REFTETRABIS15MAXN_VpF   3
#define REFTETRABIS15MAXN_FpV   7
#define REFTETRABIS15MAXN_EpF   3
#define REFTETRABIS15MAXN_FpE   4
#define REFTETRABIS15MAXN_nVpoF 5
#define REFTETRABIS15MAXN_oVpoF 3
#define REFTETRABIS15MAXN_iVpE  1
#define REFTETRABIS15MAXN_iEpF  2
#define REFTETRABIS15MAXN_nEpoE 2
#define REFTETRABIS15MAXN_nVpoE 3
#define REFTETRABIS15MAXN_nEpoF 5
#define REFTETRABIS15MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis15Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis15Desc(TShapeDesc *shape);

    // Methods
};

#endif


