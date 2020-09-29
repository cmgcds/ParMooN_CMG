// =======================================================================
// @(#)RefLineDesc.C        1.2 09/13/99
//
// Class:       TRefLineDesc
// Purpose:     refinement descriptor for line
//
// Author:      Volker Behns  07.08.97
//
// =======================================================================

#include <RefLineDesc.h>

static const Shapes DatChildType[] = { S_Line, S_Line};

static const int DatChildVertex[][LINEMAXN_VpC] = { {0, 2},  {2, 1}};
static const int DatVertexChild[][LINEMAXN_CpV] = { {0},  {1},  {0, 1}};
static const int DatVertexChildIndex[][LINEMAXN_CpV] = { {0},  {1},  {1, 0}};

static const int DatNewVertexEqOldVertex[] = { 0, 1};
static const int DatNewVertexEqOldVertexIndex[] = { 0, 1};

static const int DatInteriorVertexOfCell[] = { 2};
static const double DatPositionOfIntVert[][LINEMAXN_V] = { {0.5, 0.5}};

// Constructor
TRefLineDesc::TRefLineDesc(TShapeDesc *shape) : TRefDesc(shape)
{
  Type = LineReg;

  N_Edges = 2;
  N_Vertices = 3;
  N_Children = 2;
  N_NewVertEqOldVert = 2;


  MaxN_VpC = LINEMAXN_VpC;
  MaxN_CpV = LINEMAXN_CpV;

  ChildType = (const Shapes *) DatChildType;

  ChildVertex = (const int *) DatChildVertex;
  VertexChild = (const int *) DatVertexChild;
  VertexChildIndex = (const int *) DatVertexChildIndex;

  NewVertexEqOldVertex = (const int *) DatNewVertexEqOldVertex;
  NewVertexEqOldVertexIndex = (const int *) DatNewVertexEqOldVertexIndex;

  InteriorVertexOfCell = (const int *) DatInteriorVertexOfCell;
  PositionOfIntVert = (const double *) DatPositionOfIntVert;
}

// Methods
