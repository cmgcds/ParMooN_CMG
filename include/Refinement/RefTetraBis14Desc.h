#ifndef __REFTETRABIS14DESC__
#define __REFTETRABIS14DESC__

#include <RefDesc.h>

#define REFTETRABIS14MAXN_VpC   4
#define REFTETRABIS14MAXN_CpV   3
#define REFTETRABIS14MAXN_EpC   6
#define REFTETRABIS14MAXN_CpE   3
#define REFTETRABIS14MAXN_FpC   4
#define REFTETRABIS14MAXN_CpF   2
#define REFTETRABIS14MAXN_EpV   5
#define REFTETRABIS14MAXN_VpF   3
#define REFTETRABIS14MAXN_FpV   7
#define REFTETRABIS14MAXN_EpF   3
#define REFTETRABIS14MAXN_FpE   4
#define REFTETRABIS14MAXN_nVpoF 5
#define REFTETRABIS14MAXN_oVpoF 3
#define REFTETRABIS14MAXN_iVpE  1
#define REFTETRABIS14MAXN_iEpF  2
#define REFTETRABIS14MAXN_nEpoE 2
#define REFTETRABIS14MAXN_nVpoE 3
#define REFTETRABIS14MAXN_nEpoF 5
#define REFTETRABIS14MAXN_nFpoF 3

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraBis14Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraBis14Desc(TShapeDesc *shape);

    // Methods
};

#endif

