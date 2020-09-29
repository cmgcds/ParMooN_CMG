// =======================================================================
// @(#)IsoInterfaceJoint3D.C        1.2 09/13/99
// 
// Class:       TIsoInterfaceJoint3D
// Purpose:     connects two cells on an interface with additional
//              vertices for isoparametric reference transformation
//
// Author:      Gunar Matthies  06.08.99
//              Gunar Matthies  05.04.02
//
// =======================================================================

#include <IsoInterfaceJoint3D.h>
#include <BoundComp2D.h>
#include <BaseCell.h>

// Constructors
/** initialize the joint with the boundary parameters and one neighbour */
TIsoInterfaceJoint3D::TIsoInterfaceJoint3D(TBoundComp3D *bdcomp,
                double *param1, double *param2, TBaseCell *neigh0)
   : TInterfaceJoint3D(bdcomp, param1, param2, neigh0)
{
  ID = IsoInterfaceJoint3D;

  N_Vertices = 0;
  Vertices = NULL;
}

/** initialize the joint with the boundary parameters and two neighbours */
TIsoInterfaceJoint3D::TIsoInterfaceJoint3D(TBoundComp3D *bdcomp,
                double *param1, double *param2,
                TBaseCell *neigh0, TBaseCell *neigh1)
   : TInterfaceJoint3D(bdcomp, param1, param2, neigh0, neigh1)
{
  ID = IsoInterfaceJoint3D;

  N_Vertices = 0;
  Vertices = NULL;
}

/** initialize the joint with the boundary parameters and two neighbours */
TIsoInterfaceJoint3D::TIsoInterfaceJoint3D(TBoundComp3D *bdcomp,
                TBaseCell *neighb0)
   : TInterfaceJoint3D(bdcomp, neighb0)
{
  ID = IsoInterfaceJoint3D;

  N_Vertices = 0;
  Vertices = NULL;
}


// Destructor
TIsoInterfaceJoint3D::~TIsoInterfaceJoint3D()
   {
    if(Neighb0)
     { Neighb0 = NULL;}
  
    if(Neighb1)
     { Neighb1 = NULL;}   
   }
   

// Methods
TJoint *TIsoInterfaceJoint3D::NewInst(double newtT_0, double newT_1, TBaseCell *Me)
{
  return new TIsoInterfaceJoint3D(BoundComp, Me);
}

TJoint *TIsoInterfaceJoint3D::NewInst()
{
  return new TIsoInterfaceJoint3D(BoundComp, NULL);
}

void TIsoInterfaceJoint3D::SetVertices(int n_vertices, TVertex **vertices)
{
  delete Vertices;

  N_Vertices = n_vertices;
  Vertices = vertices;
}

void TIsoInterfaceJoint3D::GenerateVertices(int n_vertices)
{
}

