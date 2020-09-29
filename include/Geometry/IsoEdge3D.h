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

#ifndef __ISOEDGE3D__
#define __ISOEDGE3D__

#include <Vertex.h>
#include <Edge.h>


class TIsoEdge3D : public TEdge
{
  protected:
    int N_Vertices;

    TVertex **Vertices;

    int *LocDof;

    int LocVert[2];

  public:
    /** constructor with neighbours */
    TIsoEdge3D(int n_Neibs, TBaseCell **neighbs);

    /** methods */
    void SetIsoVertInfo(int n_vertices, TVertex **vertices, int *locdof, int *locvert);

    ~TIsoEdge3D();

    int GetN_Vertices()
    { return N_Vertices; }

    TVertex *GetVertex(int i)
    { return Vertices[i]; }

    int GetLocalDof (int i)
    { return LocDof[i]; }

    void GetLocalVertices(int &a, int &b)
    { a = LocVert[0]; b = LocVert[1]; }
};

#endif
