// =======================================================================
// @(#)RefHexaRegDesc.h        1.3 10/18/99
//
// Class:       TRefHexaRegDesc
// Purpose:     refinement descriptor for regular refinement of a
//              hexahedron
//
// Author:      Volker Behns  18.07.97
//
// =======================================================================

#ifndef __REFHEXAREGDESC__
#define __REFHEXAREGDESC__

#include <RefDesc.h>

#define REFHEXAREGMAXN_EpV      6
#define REFHEXAREGMAXN_FpV      12
#define REFHEXAREGMAXN_VpF      4
#define REFHEXAREGMAXN_VpC      8
#define REFHEXAREGMAXN_CpV      8
#define REFHEXAREGMAXN_EpC      12
#define REFHEXAREGMAXN_CpE      4
#define REFHEXAREGMAXN_FpC      6
#define REFHEXAREGMAXN_CpF      2
#define REFHEXAREGMAXN_EpF      4
#define REFHEXAREGMAXN_FpE      4
#define REFHEXAREGMAXN_iVpE     1
#define REFHEXAREGMAXN_oVpoF    4
#define REFHEXAREGMAXN_nVpoE    3
#define REFHEXAREGMAXN_nEpoE    2
#define REFHEXAREGMAXN_nFpoF    4
#define REFHEXAREGMAXN_nVpoF    9
#define REFHEXAREGMAXN_nEpoF    12
#define REFHEXAREGMAXN_iVpF     1
#define REFHEXAREGMAXN_iEpF     4
#define HEXAN_V                 8
#define HEXAN_E                 12
#define HEXAN_F                 6



/** refinement descriptor for regular refinement of a hexahedron */
class TRefHexaRegDesc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a hexahedron */
    TRefHexaRegDesc(TShapeDesc *shape);

    // Methods
};

#endif
