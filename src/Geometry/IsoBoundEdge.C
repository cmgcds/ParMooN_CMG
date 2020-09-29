// =======================================================================
// @(#)IsoBoundEdge.C        1.2 09/13/99
// 
// Class:       TIsoBoundEdge
// Purpose:     edge on a boundary component with additional vertices
//              for isoparametric reference transformation
//
// Author:      Gunar Matthies  06.08.99
//
// =======================================================================

#include <BoundComp2D.h>
#include <IsoBoundEdge.h>
#include <Vertex.h>

// Constructors
TIsoBoundEdge::TIsoBoundEdge(TBoundComp2D *bdcomp, double t_0, double t_1)
  : TBoundEdge(bdcomp, t_0, t_1)
{
  ID = IsoBoundEdge;

  N_Vertices = 0;
  Vertices = NULL;
}

// Methods
int TIsoBoundEdge::CheckMatchingRef(TBaseCell *Me, int J_i,
                  struct StoreGeom &Tmp)
{
  Tmp.Filled = FALSE;
  return 0;
}

#ifdef __MORTAR__

int TIsoBoundEdge::CheckMatchingRef(TBaseCell *Me, int J_i,
                  StoreGeomMortar &Tmp)
{
  Tmp.Filled = FALSE;
  return 0;
}

#endif

// create a new instance of this class
TJoint *TIsoBoundEdge::NewInst(double newT_0, double newT_1, TBaseCell *Me)
{
  return new TIsoBoundEdge(BoundComp, T_0 + newT_0*(T_1 - T_0),
                        T_0 + newT_1*(T_1 - T_0));
}

TJoint *TIsoBoundEdge::NewInst()
{
  return new TIsoBoundEdge(BoundComp, T_0, T_1);
}

void TIsoBoundEdge::SetVertices(int n_vertices, TVertex **vertices)
{
  if(Vertices)
    delete Vertices;

  N_Vertices = n_vertices;
  Vertices = vertices;
}

void TIsoBoundEdge::GenerateVertices(int n_vertices)
{
  if(N_Vertices != n_vertices)
  {
    if(Vertices)
      delete Vertices;

    N_Vertices = n_vertices;
    Vertices = new TVertex*[N_Vertices];
  }
#ifdef __2D__
  double t, x, y;
  for(int i=0;i<N_Vertices;i++)
  {
    t = T_0 + ((i+1)*(T_1-T_0))/(N_Vertices+1);
    BoundComp->GetXYofT(t, x, y);
    Vertices[i] = new TVertex(x, y);
  } // endfor i
#endif // __2D__
}

void TIsoBoundEdge::GeneratemidVert(int n_vertices, double*X, double*Y)
{
   if(Vertices)
     delete Vertices;

  N_Vertices = n_vertices;
  Vertices = new TVertex*[N_Vertices];

#ifdef __2D__
  double  x, y;
  for(int i=0;i<N_Vertices;i++)
  {
    x = (X[0]+X[1])/2.0;
    y = (Y[0]+Y[1])/2.0;
    
    //double r = sqrt(x*x + y*y);
    
    //cout << r<< " x " << x << " y " << y <<endl;
    //x /=r;
    //y /=r;
    //cout << " x " << x << " y " << y <<endl;
    Vertices[i] = new TVertex(x, y);;
  } // endfor i
#endif // __2D__
}

void TIsoBoundEdge::DeleteVertices()
{
  if(Vertices)
    delete Vertices;
  N_Vertices = 0;
  Vertices = NULL;
 }
 
 
