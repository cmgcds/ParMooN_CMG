// ***********************************************************************
// Q0 element, discontinuous, 3D
// ***********************************************************************

// number of degrees of freedom
static int C_H_Q0_3D_NDOF = 1;

// number of dofs on the closure of the joints
static int C_H_Q0_3D_JointDOF = 0;

// which local dofs are on the joints
static int *C_H_Q0_3D_J0 = NULL;
static int *C_H_Q0_3D_J1 = NULL;
static int *C_H_Q0_3D_J2 = NULL;
static int *C_H_Q0_3D_J3 = NULL;
static int *C_H_Q0_3D_J4 = NULL;
static int *C_H_Q0_3D_J5 = NULL;

static int *C_H_Q0_3D_J[6] = { C_H_Q0_3D_J0, C_H_Q0_3D_J1,
                               C_H_Q0_3D_J2, C_H_Q0_3D_J3,
                               C_H_Q0_3D_J4, C_H_Q0_3D_J5};

// number of inner dofs
static int C_H_Q0_3D_NInner = 1;

// array containing the numbers for the inner dofs
static int C_H_Q0_3D_Inner[1] = { 0 };

#ifdef _MPI   
// number of dofs on the closure of the edges
static int C_H_Q0_3D_EdgeDOF = 0;

// which local dofs are on the joints
static int *C_H_Q0_3D_E0 =  NULL;
static int *C_H_Q0_3D_E1 =  NULL;
static int *C_H_Q0_3D_E2 =  NULL;
static int *C_H_Q0_3D_E3 =  NULL;

static int *C_H_Q0_3D_E4 =  NULL;
static int *C_H_Q0_3D_E5 =  NULL;
static int *C_H_Q0_3D_E6 =  NULL;
static int *C_H_Q0_3D_E7 =  NULL;

static int *C_H_Q0_3D_E8 =  NULL;
static int *C_H_Q0_3D_E9 =  NULL;
static int *C_H_Q0_3D_E10 =  NULL;
static int *C_H_Q0_3D_E11 =  NULL;

static int *C_H_Q0_3D_E[12] = { C_H_Q0_3D_E0, C_H_Q0_3D_E1, C_H_Q0_3D_E2, C_H_Q0_3D_E3,
                                C_H_Q0_3D_E4, C_H_Q0_3D_E5, C_H_Q0_3D_E6, C_H_Q0_3D_E7,
                                C_H_Q0_3D_E8, C_H_Q0_3D_E9, C_H_Q0_3D_E10, C_H_Q0_3D_E11};

// number of dofs on the closure of the vertices
static int C_H_Q0_3D_VertDOF = 0;

// array containing the numbers for the vertices dofs
static int *C_H_Q0_3D_Vert = NULL;

#endif


static char C_H_Q0_3D_String[] = "C_H_Q0_3D";

TFEDesc3D *FE_C_H_Q0_3D_Obj=new TFEDesc3D(C_H_Q0_3D_String, C_H_Q0_3D_NDOF, 
                                C_H_Q0_3D_JointDOF,
                                C_H_Q0_3D_J, C_H_Q0_3D_NInner, C_H_Q0_3D_Inner
#ifdef _MPI
                                ,C_H_Q0_3D_EdgeDOF,  C_H_Q0_3D_E, C_H_Q0_3D_VertDOF,
                                 C_H_Q0_3D_Vert
#endif
                                 );

