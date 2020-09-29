// =======================================================================
// @(#)RefQuadBis1Desc.h        1.1 02/08/99
//
// Class:       TRefQuadBis1Desc
// Purpose:     refinement descriptor for a bisection of a quadrangle
//              based on edge 1
//
// Author:      Volker Behns  08.02.99
//
// =======================================================================

#ifndef __REFQUADBIS1DESC__
#define __REFQUADBIS1DESC__

#include <RefDesc.h>

#ifndef __REFQUADBIS0DESC__
  #define QUADBIN_E         4
  #define QUADBIMAXN_VpC    4
  #define QUADBIMAXN_CpV    2
  #define QUADBIMAXN_EpC    4
  #define QUADBIMAXN_CpE    2
  #define QUADBIMAXN_EpV    3
  #define QUADBIMAXN_iVpE   1
  #define QUADBIMAXN_nVpoE  3
  #define QUADBIMAXN_nEpoE  2
#endif

/** refinement descriptor for a bisection of a quadrangle based on edge 1 */
class TRefQuadBis1Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for bisection of a quadrangle */
    TRefQuadBis1Desc(TShapeDesc *shape);

    // Methods
};

#endif
