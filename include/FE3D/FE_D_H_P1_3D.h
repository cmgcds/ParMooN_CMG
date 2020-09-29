// ***********************************************************************
// P1 element, discontinuous, 3D
// ***********************************************************************

// number of degrees of freedom
static int D_H_P1_3D_NDOF = 4;

// number of dofs on the closure of the joints
static int D_H_P1_3D_JointDOF = 0;

// which local dofs are on the joints
static int *D_H_P1_3D_J0 = NULL;
static int *D_H_P1_3D_J1 = NULL;
static int *D_H_P1_3D_J2 = NULL;
static int *D_H_P1_3D_J3 = NULL;
static int *D_H_P1_3D_J4 = NULL;
static int *D_H_P1_3D_J5 = NULL;

static int *D_H_P1_3D_J[6] = { D_H_P1_3D_J0, D_H_P1_3D_J1,
                             D_H_P1_3D_J2, D_H_P1_3D_J3,
                             D_H_P1_3D_J4, D_H_P1_3D_J5};

#ifdef _MPI   
// number of dofs on the closure of the edges
static int D_H_P1_3D_EdgeDOF = 0;

static int *D_H_P1_3D_E0 = NULL;
static int *D_H_P1_3D_E1 = NULL;
static int *D_H_P1_3D_E2 = NULL;
static int *D_H_P1_3D_E3 = NULL;
static int *D_H_P1_3D_E4 = NULL;
static int *D_H_P1_3D_E5 = NULL;
static int *D_H_P1_3D_E6 = NULL;
static int *D_H_P1_3D_E7 = NULL;
static int *D_H_P1_3D_E8 = NULL;
static int *D_H_P1_3D_E9 = NULL;
static int *D_H_P1_3D_E10 = NULL;
static int *D_H_P1_3D_E11 = NULL;

static int *D_H_P1_3D_E[12] = {D_H_P1_3D_E0, D_H_P1_3D_E1, D_H_P1_3D_E2, D_H_P1_3D_E3,
                               D_H_P1_3D_E4, D_H_P1_3D_E5, D_H_P1_3D_E6, D_H_P1_3D_E7,
                               D_H_P1_3D_E8, D_H_P1_3D_E9, D_H_P1_3D_E10, D_H_P1_3D_E11  };


// number of dofs on the closure of the vertex
static int D_H_P1_3D_VertDOF = 0;

// array containing the numbers for the vertices dofs
static int *D_H_P1_3D_Vert = NULL;

#endif


// number of inner dofs
static int D_H_P1_3D_NInner = 4;

// array containing the numbers for the inner dofs (here is no inner dof)
static int D_H_P1_3D_Inner[] = { 0, 1, 2, 3 };

static char D_H_P1_3D_String[] = "D_H_P1_3D";

TFEDesc3D *FE_D_H_P1_3D_Obj=new TFEDesc3D(D_H_P1_3D_String, D_H_P1_3D_NDOF, 
                                D_H_P1_3D_JointDOF,
                                D_H_P1_3D_J, D_H_P1_3D_NInner,
                                D_H_P1_3D_Inner
#ifdef _MPI
                                , D_H_P1_3D_EdgeDOF, D_H_P1_3D_E,
                                D_H_P1_3D_VertDOF, D_H_P1_3D_Vert
#endif
                                 );
