// =======================================================================
// @(#)RefMortar0Desc.h        1.1 10/30/98
//
// Class:       TRefMortar0Desc
// Purpose:     refinement descriptor for a mortar cell (i.e. generate a
//              Nx1 grid with base edge 0)
//
// Author:      Volker Behns  30.12.97
//
// =======================================================================

#ifndef __REFMORTAR0DESC__
#define __REFMORTAR0DESC__

#include <RefDesc.h>

#ifndef __REFMORTAR1DESC__
  #define MORTARRMAXN_VpC    4
  #define MORTARRMAXN_CpV    2
  #define MORTARRMAXN_EpC    4
  #define MORTARRMAXN_CpE    2
  #define MORTARRMAXN_EpV    3
#endif

/** refinement descriptor for a mortar cell */
class TRefMortar0Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for mortar refinement of a quadrangle */
    TRefMortar0Desc(TShapeDesc *shape, int Mortar_Ni, int N);

    // Methods
};

#endif
