// ***********************************************************************
// P1 element, discontinuous, 3D
// 
// Author:     Sashikumaar Ganesan
//
// ***********************************************************************

// number of degrees of freedom
static int D_T_P1_3D_NDOF = 4;

// number of dofs on the closure of the joints
static int D_T_P1_3D_JointDOF = 0;

// which local dofs are on the joints
static int *D_T_P1_3D_J0 = NULL;
static int *D_T_P1_3D_J1 = NULL;
static int *D_T_P1_3D_J2 = NULL;
static int *D_T_P1_3D_J3 = NULL;

static int *D_T_P1_3D_J[4] = { D_T_P1_3D_J0, D_T_P1_3D_J1,
                               D_T_P1_3D_J2, D_T_P1_3D_J3 };

// number of inner dofs
#ifdef _MPI   
// number of dofs on the closure of the edges
static int D_T_P1_3D_EdgeDOF = 0;

static int *D_T_P1_3D_E0 = NULL;
static int *D_T_P1_3D_E1 = NULL;
static int *D_T_P1_3D_E2 = NULL;
static int *D_T_P1_3D_E3 = NULL;
static int *D_T_P1_3D_E4 = NULL;
static int *D_T_P1_3D_E5 = NULL;


static int *D_T_P1_3D_E[6] = {D_T_P1_3D_E0, D_T_P1_3D_E1, D_T_P1_3D_E2, D_T_P1_3D_E3,
                               D_T_P1_3D_E4, D_T_P1_3D_E5 };


// number of dofs on the closure of the vertex
static int D_T_P1_3D_VertDOF = 0;

// array containing the numbers for the vertices dofs
static int *D_T_P1_3D_Vert = NULL;

#endif
			       
			       
static int D_T_P1_3D_NInner = 4;

// array containing the numbers for the inner dofs 
static int D_T_P1_3D_Inner[] = { 0, 1, 2, 3 };

static char D_T_P1_3D_String[] = "D_T_P1_3D";

TFEDesc3D *FE_D_T_P1_3D_Obj=new TFEDesc3D(D_T_P1_3D_String, D_T_P1_3D_NDOF,
                                D_T_P1_3D_JointDOF,
                                D_T_P1_3D_J, D_T_P1_3D_NInner,
                                D_T_P1_3D_Inner
#ifdef _MPI
                                , D_T_P1_3D_EdgeDOF, D_T_P1_3D_E,
                                D_T_P1_3D_VertDOF, D_T_P1_3D_Vert
#endif
                                 );