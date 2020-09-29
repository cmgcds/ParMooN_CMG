// ***********************************************************************
// P2 element, nonconforming, 3D
// ***********************************************************************

// number of degrees of freedom
static int N_T_P2_3D_NDOF = 13;

// number of dofs on the closure of the joints
static int N_T_P2_3D_JointDOF = 3;

// which local dofs are on the joints
static int N_T_P2_3D_J0[3] = { 0,  1,  2 };
static int N_T_P2_3D_J1[3] = { 3,  4,  5 };
static int N_T_P2_3D_J2[3] = { 6,  7,  8 };
static int N_T_P2_3D_J3[3] = { 9, 10, 11 };

static int *N_T_P2_3D_J[4] = { N_T_P2_3D_J0, N_T_P2_3D_J1,
                             N_T_P2_3D_J2, N_T_P2_3D_J3 };

// number of inner dofs
static int N_T_P2_3D_NInner = 1;

// array containing the numbers for the inner dofs (here is no inner dof)
static int N_T_P2_3D_Inner[1] = { 12 };

static char N_T_P2_3D_String[] = "N_T_P2_3D";

TFEDesc3D *FE_N_T_P2_3D_Obj=new TFEDesc3D(N_T_P2_3D_String, N_T_P2_3D_NDOF, 
                                N_T_P2_3D_JointDOF,
                                N_T_P2_3D_J, N_T_P2_3D_NInner, N_T_P2_3D_Inner);
