// =======================================================================
// @(#)Line.C        1.1 10/30/98
//
// Class:       TLine
// Purpose:     shape descriptor of a line
//
// Author:      Volker Behns  07.08.97
//
// =======================================================================

#include <Line.h>

static const int DatEdgeVertex[][2] = { {0, 1}};
static const int DatVertexEdge[][LINEMAXN_EpV] = { {0},  {0}};

// Constructor
TLine::TLine()
{
  MaxN_EpV = LINEMAXN_EpV;

  EdgeVertex = (const int *) DatEdgeVertex;
  VertexEdge = (const int *) DatVertexEdge;

  Type = S_Line;
  N_Vertices = 2;
  N_Edges = 2;
  N_Joints = 2;
}

// Methods
double TLine::GetMeasure(TVertex **Verts)
{
  double x1,x2,y1,y2;

  x1 = Verts[0]->GetX();
  y1 = Verts[0]->GetY();
  x2 = Verts[1]->GetX();
  y2 = Verts[1]->GetY();

  return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}
