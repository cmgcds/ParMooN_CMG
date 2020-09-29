// =======================================================================
// @(#)RefQuad1Conf1Desc.C        1.2 09/13/99
//
// Class:       TRefQuad1Conf1Desc
// Purpose:     refinement descriptor for conforming closure of a
//              quadrangle with 1 hanging node
//
// Author:      Matthias Ebeling  26.08.99
//
// =======================================================================

#include <RefQuad1Conf1Desc.h>

static const Shapes DatChildType[] = { Triangle, Triangle,
                                       Triangle };

static const Refinements DatEdgeType[] = { NoRef, NoRef, NoRef, LineReg};

static const int DatChildVertex[][QUAD1ConfMAXN_VpC] =
                 { {0, 1, 4},  {1, 2, 4},  {0, 4, 3}};
static const int DatVertexChild[][QUAD1ConfMAXN_CpV] =
                 { {0, 2},  {0, 1},  {1},  {2},  {0, 1, 2}};
                 
static const int DatVertexChildIndex[][QUAD1ConfMAXN_CpV] =
               { {0, 0},  {1, 0},  {1},  {2},  {2, 2, 1}};

static const int DatVertexChildLen[] =
               { 2,  2,  1,  1,  3};

static const int DatChildEdge[][QUAD1ConfMAXN_EpC] =
               { {0, 6, 5},  {1, 2, 6},  {5, 3, 4}};
static const int DatEdgeChild[][QUAD1ConfMAXN_CpE] =
               { {0},  {1},  {1},  {2},  {2},  {0, 2},  {0, 1}};
static const int DatEdgeChildIndex[][QUAD1ConfMAXN_CpE] = 
               { {0},  {0},  {1},  {1},  {2},  {2, 0},  {1, 2}};
static const int DatEdgeChildLen[] =
               { 1,  1,  1,  1,  1,  2,  2};

static const int DatEdgeVertex[][2] =
               { {0, 1},  {1, 2},  {2, 4},  {4, 3},  {3, 0},  {0, 4},
                 {4, 1}};
static const int DatVertexEdge[][QUAD1ConfMAXN_EpV] = 
               { {0, 5, 4},  {0, 1, 6},  {1, 2},  {3, 4}, 
                 {2, 3, 5, 6}};
static const int DatVertexEdgeIndex[][QUAD1ConfMAXN_EpV] =
               { {0, 0, 1},  {1, 0, 1},  {1, 0},  {1, 0},
                 {1, 0, 1, 0}};
static const int DatVertexEdgeLen[] =
               { 3,  3,  2,  2,  4};

static const int DatNewVertexEqOldVertex[] = { 0, 1, 2, 3};
static const int DatNewVertexEqOldVertexIndex[] = { 1, 2, 3, 0};

                                 
static const int DatInteriorEdgeOfCell[] = { 5, 6};
static const int DatInteriorVertexOfEdge[][QUAD1ConfMAXN_iVpE] =
               { {-1}, {-1}, {-1}, {4}};
static const int DatInteriorVertexOfEdgeLen[] = { 0, 0, 0, 1};

static const int DatOldEdgeNewVertex[][QUAD1ConfMAXN_nVpoE] =
               { {3, 0},  {0, 1},  {1, 2},  {2, 4, 3}};

static const int DatOldEdgeNewEdge[][QUAD1ConfMAXN_nEpoE] =
               { {4},  {0},  {1},  {2, 3}};

static const int DatOldEdgeNewLocEdge[][QUADConfN_E] =
               { {-1, 0, -1, -1}, {-1, -1, 0, 1},
                 {2, -1, -1, 1} };

static const int DatNewEdgeOldEdge[] =
               { 1,  2,  3,  3,  0, -1,  -1};

static const int DatNewEdgeEqOldEdge[] = { 0, 1, 4};
static const int DatNewEdgeEqOldEdgeIndex[] = { 1, 2, 0};

// Constructor
TRefQuad1Conf1Desc::TRefQuad1Conf1Desc(TShapeDesc *shape) : TRefDesc(shape)
{
  Type = Quad1Conf1;

  // set all numbers
  N_Edges = 7;
  N_Vertices = 5;
  N_Children = 3;
  N_NewVertEqOldVert = 4;
  N_NewEdgeEqOldEdge = 3;
  N_InnerVertices = 0;    
  N_InnerEdges = 2;

  // initialize all dimension values
  MaxN_VpC = QUAD1ConfMAXN_VpC;
  MaxN_CpV = QUAD1ConfMAXN_CpV;
  MaxN_EpC = QUAD1ConfMAXN_EpC;
  MaxN_CpE = QUAD1ConfMAXN_CpE;
  MaxN_EpV = QUAD1ConfMAXN_EpV;
  MaxN_iVpE = QUAD1ConfMAXN_iVpE;
  MaxN_nVpoE = QUAD1ConfMAXN_nVpoE;
  MaxN_nEpoE = QUAD1ConfMAXN_nEpoE;

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


  InteriorEdgeOfCell = (const int *) DatInteriorEdgeOfCell;
  InteriorVertexOfEdge = (const int *) DatInteriorVertexOfEdge;
  InteriorVertexOfEdgeLen = (const int *) DatInteriorVertexOfEdgeLen;

  OldEdgeNewVertex = (const int *) DatOldEdgeNewVertex;

  OldEdgeNewEdge = (const int *) DatOldEdgeNewEdge;
  OldEdgeNewLocEdge = (const int *) DatOldEdgeNewLocEdge;
  NewEdgeOldEdge = (const int *) DatNewEdgeOldEdge;

  NewEdgeEqOldEdge = (const int *) DatNewEdgeEqOldEdge;
  NewEdgeEqOldEdgeIndex = (const int *) DatNewEdgeEqOldEdgeIndex;
}

// Methods
