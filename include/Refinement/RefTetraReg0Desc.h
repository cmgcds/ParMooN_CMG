// =======================================================================
// @(#)RefTetraReg0Desc.h        1.2 10/18/99
//
// Class:       TRefTetraRegDesc
// Purpose:     refinement descriptor for regular refinement of a
//              tetrahedron (variant 0)
//
// Authors:      Volker Behns  18.07.97
//               Matthias Ebelibg 06.09.99
//
// =======================================================================

#ifndef __REFTETRAREG0DESC__
#define __REFTETRAREG0DESC__

#include <RefDesc.h>

#define REFTETRAREGMAXN_EpV      7
#define REFTETRAREGMAXN_FpV     12
#define REFTETRAREGMAXN_VpF      3
#define REFTETRAREGMAXN_VpC      4
#define REFTETRAREGMAXN_CpV      6
#define REFTETRAREGMAXN_EpC      6
#define REFTETRAREGMAXN_CpE      4
#define REFTETRAREGMAXN_FpC      4
#define REFTETRAREGMAXN_CpF      2
#define REFTETRAREGMAXN_EpF      3
#define REFTETRAREGMAXN_FpE      4
#define REFTETRAREGMAXN_iVpE     1
#define REFTETRAREGMAXN_nVpoE    3
#define REFTETRAREGMAXN_nEpoE    2
#define REFTETRAREGMAXN_nFpoF    4
#define REFTETRAREGMAXN_nVpoF    6
#define REFTETRAREGMAXN_oVpoF    3
#define REFTETRAREGMAXN_nEpoF    9
#define REFTETRAREGMAXN_iEpF     3
#define TETRAN_V                 4
#define TETRAN_E                 6
#define TETRAN_F                 4

/** refinement descriptor for regular refinement of a tetrahedron */
class TRefTetraReg0Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a tetrahedron */
    TRefTetraReg0Desc(TShapeDesc *shape);

    // Methods
};

#endif
