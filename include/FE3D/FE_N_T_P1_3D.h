// ***********************************************************************
// P1 element, nonconforming, 3D
// ***********************************************************************

// number of degrees of freedom
static int N_T_P1_3D_NDOF = 4;

// number of dofs on the closure of the joints
static int N_T_P1_3D_JointDOF = 1;

// which local dofs are on the joints
static int N_T_P1_3D_J0[1] = { 0 };
static int N_T_P1_3D_J1[1] = { 1 };
static int N_T_P1_3D_J2[1] = { 2 };
static int N_T_P1_3D_J3[1] = { 3 };

static int *N_T_P1_3D_J[4] = { N_T_P1_3D_J0, N_T_P1_3D_J1,
                             N_T_P1_3D_J2, N_T_P1_3D_J3 };

// number of inner dofs
static int N_T_P1_3D_NInner = 0;

// array containing the numbers for the inner dofs (here is no inner dof)
static int *N_T_P1_3D_Inner = NULL;

#ifdef _MPI   
// number of dofs on the closure of the edges
static int N_T_P1_3D_EdgeDOF = 0;

// which local dofs are on the joints
static int *N_T_P1_3D_E0 = NULL;
static int *N_T_P1_3D_E1 = NULL;
static int *N_T_P1_3D_E2 = NULL;
static int *N_T_P1_3D_E3 = NULL;
static int *N_T_P1_3D_E4 = NULL;
static int *N_T_P1_3D_E5 = NULL;


static int *N_T_P1_3D_E[6] = { N_T_P1_3D_E0, N_T_P1_3D_E1, N_T_P1_3D_E2, N_T_P1_3D_E3,
                               N_T_P1_3D_E4, N_T_P1_3D_E5};

// number of dofs on the closure of the vertices
static int N_T_P1_3D_VertDOF = 0;

// array containing the numbers for the vertices dofs
static int *N_T_P1_3D_Vert = NULL;

#endif

static char N_T_P1_3D_String[] = "N_T_P1_3D";

TFEDesc3D *FE_N_T_P1_3D_Obj=new TFEDesc3D(N_T_P1_3D_String, N_T_P1_3D_NDOF, 
                                N_T_P1_3D_JointDOF,
                                N_T_P1_3D_J, N_T_P1_3D_NInner, N_T_P1_3D_Inner
#ifdef _MPI
                                ,N_T_P1_3D_EdgeDOF,  N_T_P1_3D_E, N_T_P1_3D_VertDOF,
                                 N_T_P1_3D_Vert
#endif 
);
