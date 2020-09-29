// =======================================================================
// @(#)RefTriBis0Desc.C        1.3 09/13/99
//
// Class:       TRefTriBis0Desc
// Purpose:     bisection of edge 0
//              into two triangles
//
// Author:      Volker Behns  16.03.98
//              Gunar Matthies 03.09.1998
//
// =======================================================================

#include <RefTriBis21Desc.h>

static const Shapes DatChildType[] = { Triangle, Triangle, Triangle};

static const Refinements DatEdgeType[] = { NoRef, LineReg, LineReg };

static const int DatChildVertex[][TRIBI21MAXN_VpC] = {{0,1,3},{2,3,4},{1,4,3}};
static const int DatVertexChild[][TRIBI21MAXN_CpV] = {{0},{0,2},{1},{0,1,2},{1,2}};
static const int DatVertexChildIndex[][TRIBI21MAXN_CpV] = {{0},{1,0},{0},{2,1,2},{2,1}};
static const int DatVertexChildLen[] = {1,2,1,3,2};

static const int DatChildEdge[][TRIBI21MAXN_EpC] = {{0,5,4},{3,6,2},{1,6,5}};
static const int DatEdgeChild[][TRIBI21MAXN_CpE] = {{0},{2},{1},{1},{0},{0,2},{1,2}};
static const int DatEdgeChildIndex[][TRIBI21MAXN_CpE] = {{0},{0},{2},{0},{2},{1,2},{1,1}};
static const int DatEdgeChildLen[] = {1,1,1,1,1,2,2};

static const int DatEdgeVertex[][2] = {{0,1},{1,4},{4,2},{2,3},{3,0},{3,1},{3,4}};
static const int DatVertexEdge[][TRIBI21MAXN_EpV] = {{0,4},{0,1,5},{2,3},{3,4,5,6},{1,2,6}};
static const int DatVertexEdgeIndex[][TRIBI21MAXN_EpV] = {{0,1},{1,0,1},{1,0},{1,0,0,0},{1,0,1}};
static const int DatVertexEdgeLen[] = {2,3,2,4,3};

static const int DatNewVertexEqOldVertex[] = { 0, 1, 2};
static const int DatNewVertexEqOldVertexIndex[] = { 0, 1, 2};

static const int DatNewEdgeEqOldEdge[] = {0};
static const int DatNewEdgeEqOldEdgeIndex[] = {0};

static const int DatInteriorVertexOfEdge[][TRIBI21MAXN_iVpE] = { {}, {4}, {3}};
static const int DatInteriorVertexOfEdgeLen[] = { 0, 1, 1};

static const int DatInteriorEdgeOfCell[] = {5, 6};

static const int DatOldEdgeNewVertex[][TRIBI21MAXN_nVpoE] = { {0, 1}, { 1, 4, 2}, {2, 3, 0}};

static const int DatOldEdgeNewEdge[][TRIBI21MAXN_nEpoE] =
               { {0}, {1, 2}, {3, 4}};

static const int DatOldEdgeNewLocEdge[][TRIBI21N_E] =
               { {0, -1, 2}, {-1, 2, 0}, {-1, 0, -1} };

static const int DatNewEdgeOldEdge[] =
               { 0, 1, 1, 2, 2, -1, -1};

// Constructor
TRefTriBis21Desc::TRefTriBis21Desc(TShapeDesc *shape) : TRefDesc(shape)
{
  Type = TriBis21;

  // set all numbers
  N_Edges = 7;
  N_Vertices = 5;
  N_Children = 3;
  N_NewVertEqOldVert = 3;
  N_NewEdgeEqOldEdge = 1;
  N_InnerEdges = 2;

  // initialize all dimension values
  MaxN_VpC = TRIBI21MAXN_VpC;
  MaxN_CpV = TRIBI21MAXN_CpV;
  MaxN_EpC = TRIBI21MAXN_EpC;
  MaxN_CpE = TRIBI21MAXN_CpE;
  MaxN_EpV = TRIBI21MAXN_EpV;
  MaxN_nVpoE = TRIBI21MAXN_nVpoE;
  MaxN_nEpoE = TRIBI21MAXN_nEpoE;

  // initialize all pointers
  ChildType = (const Shapes *) DatChildType;
  EdgeType = (const Refinements *) DatEdgeType;

  ChildVertex = (const int *) DatChildVertex;
  VertexChild = (const int *) DatVertexChild;
  VertexChildIndex = (const int *) DatVertexChildIndex;
  VertexChildLen = (const int *) DatVertexChildLen;

  ChildEdge = (const int *) DatChildEdge;
  EdgeChild = (const int *) DatEdgeChild;
  EdgeChildIndex = (const int *) DatEdgeChildIndex;
  EdgeChildLen = (const int *) DatEdgeChildLen;

  EdgeVertex = (const int *) DatEdgeVertex;
  VertexEdge = (const int *) DatVertexEdge;
  VertexEdgeIndex = (const int *) DatVertexEdgeIndex;
  VertexEdgeLen = (const int *) DatVertexEdgeLen;

  NewVertexEqOldVertex = (const int *) DatNewVertexEqOldVertex;
  NewVertexEqOldVertexIndex = (const int *) DatNewVertexEqOldVertexIndex;

  NewEdgeEqOldEdge = (const int *) DatNewEdgeEqOldEdge;
  NewEdgeEqOldEdgeIndex = (const int *) DatNewEdgeEqOldEdgeIndex;

  InteriorEdgeOfCell = (const int *) DatInteriorEdgeOfCell;

  OldEdgeNewVertex = (const int *) DatOldEdgeNewVertex;

  OldEdgeNewEdge = (const int *) DatOldEdgeNewEdge;
  OldEdgeNewLocEdge = (const int *) DatOldEdgeNewLocEdge;
  NewEdgeOldEdge = (const int *) DatNewEdgeOldEdge;

  InteriorVertexOfEdge = (const int *) DatInteriorVertexOfEdge;
  InteriorVertexOfEdgeLen = (const int *) DatInteriorVertexOfEdgeLen;
}

// Methods
