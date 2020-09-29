// =======================================================================
// @(#)IsoEdge3D
// 
// Class:       TIsoEdge3D
// Purpose:     class for iso edges in 3D
//
// Author:      Sashikumaar Ganesan  06.09.2010
//
// History:
//
// ======================================================================= 


#include <IsoEdge3D.h>

/** constructor with neighbours */
TIsoEdge3D::TIsoEdge3D(int n_Neibs, TBaseCell **neighbs):TEdge(n_Neibs, neighbs)
{
  EdgeID = IsoEdge3D;
}


void TIsoEdge3D::SetIsoVertInfo (int n_vertices, TVertex **vertices, int *locdof, int *locvert)
{
   N_Vertices = n_vertices;

   Vertices = new TVertex* [n_vertices];
   LocDof   = new int [n_vertices];

   for (int i=0;i<n_vertices;++i)
   {
     Vertices[i] = vertices[i];
     LocDof[i] = locdof[i];
   }

   LocVert[0] = locvert[0];
   LocVert[1] = locvert[1];
}

// DTOR
TIsoEdge3D::~TIsoEdge3D()
{
  if ( Vertices ) delete [] Vertices;
  if ( LocDof ) delete [] LocDof;
}
