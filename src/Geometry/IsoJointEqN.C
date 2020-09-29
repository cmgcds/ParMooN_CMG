// =======================================================================
// @(#)IsoJointEqN.C        1.2 10/18/99
// 
// Class:       TIsoJointEqN
// Purpose:     connects two cells
//              has additional vertices for isoparametric reference 
//              transformation
//
// Author:      Gunar Matthies 06.08.1999
//
// History:     start of implementation (Gunar Matthies 06.08.99)
//
// =======================================================================

#include <BaseCell.h>
#include <Constants.h>
#include <IsoJointEqN.h>
#include <RefDesc.h>

// Constructors
TIsoJointEqN::TIsoJointEqN(TBaseCell *neighb0) : TJointEqN(neighb0)
{
  ID = IsoJointEqN;

  N_Vertices = 0;
  Vertices = NULL;
}

TIsoJointEqN::TIsoJointEqN(TBaseCell *neighb0, TBaseCell *neighb1) 
  : TJointEqN(neighb0, neighb1)
{
  ID = IsoJointEqN;

  N_Vertices = 0;
  Vertices = NULL;
}

// Methods
void TIsoJointEqN::SetVertices(int n_vertices, TVertex **vertices)
{
  N_Vertices = n_vertices;

  Vertices = vertices;
}
