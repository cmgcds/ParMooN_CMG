// =======================================================================
// @(#)Triangle.C        1.1 10/30/98
//
// Class:       TTriangle
// Purpose:     shape descriptor of a triangle
//
// Author:      Volker Behns  30.07.97
//
// =======================================================================

#include <Triangle.h>
#include <Constants.h>
#include <math.h>

static const int DatEdgeVertex[][2] = { {0, 1},  {1, 2},  {2, 0}};
static const int DatVertexEdge[][TRIMAXN_EpV] = { {2, 0},  {0, 1},  {1, 2}};

// Constructor
TTriangle::TTriangle()
{
  MaxN_EpV = TRIMAXN_EpV;

  EdgeVertex = (const int *) DatEdgeVertex;
  VertexEdge = (const int *) DatVertexEdge;

  Type = Triangle;
  N_Vertices = 3;
  N_Edges = 3;
  N_Joints = 3;
}

// Methods
double TTriangle::GetDiameter(TVertex **Verts)
{
  double diffX1 = Verts[0]->GetX() - Verts[1]->GetX();
  double diffY1 = Verts[0]->GetY() - Verts[1]->GetY();
  double diffX2 = Verts[1]->GetX() - Verts[2]->GetX();
  double diffY2 = Verts[1]->GetY() - Verts[2]->GetY();
  double diffX3 = Verts[2]->GetX() - Verts[0]->GetX();
  double diffY3 = Verts[2]->GetY() - Verts[0]->GetY();

  return sqrt(MAX(diffX1*diffX1 + diffY1*diffY1,
              MAX(diffX2*diffX2 + diffY2*diffY2,
                  diffX3*diffX3 + diffY3*diffY3)));
}

double TTriangle::GetShortestEdge(TVertex **Verts)
{
  double diffX1 = Verts[0]->GetX() - Verts[1]->GetX();
  double diffY1 = Verts[0]->GetY() - Verts[1]->GetY();
  double diffX2 = Verts[1]->GetX() - Verts[2]->GetX();
  double diffY2 = Verts[1]->GetY() - Verts[2]->GetY();
  double diffX3 = Verts[2]->GetX() - Verts[0]->GetX();
  double diffY3 = Verts[2]->GetY() - Verts[0]->GetY();
  double len = 1e10;
  
  if (sqrt(diffX1*diffX1 + diffY1*diffY1)< len)
      len = sqrt(diffX1*diffX1 + diffY1*diffY1);
  if (sqrt(diffX2*diffX2 + diffY2*diffY2)< len)
      len = sqrt(diffX2*diffX2 + diffY2*diffY2);
  if (sqrt(diffX3*diffX3 + diffY3*diffY3)< len)
      len = sqrt(diffX3*diffX3 + diffY3*diffY3);

  return(len);
}

double TTriangle::GetLengthWithReferenceMap(TVertex **Verts)
{
  double x0, x1, x2, y0, y1, y2;
  double xc1, xc2, yc1, yc2;

  x0 = Verts[0]->GetX();
  y0 = Verts[0]->GetY();
  x1 = Verts[1]->GetX();
  y1 = Verts[1]->GetY();
  x2 = Verts[2]->GetX();
  y2 = Verts[2]->GetY();

  //double xc0=x0;
  xc1=x1-x0;
  xc2=x2-x0;

  //double yc0=y0;
  yc1=y1-y0;
  yc2=y2-y0;

  return sqrt(fabs(xc1*yc2-xc2*yc1));
}

double TTriangle::GetMeasure(TVertex **Verts)
{
  double x1,x2,x3,y1,y2,y3;

  x1 = Verts[0]->GetX();
  y1 = Verts[0]->GetY();
  x2 = Verts[1]->GetX();
  y2 = Verts[1]->GetY();
  x3 = Verts[2]->GetX();
  y3 = Verts[2]->GetY();

  return 0.5*fabs(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);
}
