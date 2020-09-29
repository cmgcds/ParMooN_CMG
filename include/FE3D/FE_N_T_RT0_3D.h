// ***********************************************************************
// Raviart-Thomas element of zero-th order on tetrahedra, 3D
// ***********************************************************************

// number of degrees of freedom
static int N_T_RT0_3D_NDOF = 4;

// number of dofs on the closure of each joints
static int N_T_RT0_3D_JointDOF = 1;

// which local dofs are on the joints
static int N_T_RT0_3D_J0[1] = { 0 };
static int N_T_RT0_3D_J1[1] = { 1 };
static int N_T_RT0_3D_J2[1] = { 2 };
static int N_T_RT0_3D_J3[1] = { 3 };

static int *N_T_RT0_3D_J[4] = { N_T_RT0_3D_J0, N_T_RT0_3D_J1,
                                N_T_RT0_3D_J2, N_T_RT0_3D_J3 };

// number of inner dofs
static int N_T_RT0_3D_NInner = 0;

// array containing the numbers for the inner dofs (here is no inner dof)
static int *N_T_RT0_3D_Inner = NULL;

// number of outer dofs
static int N_T_RT0_3D_NOuter = 4;

// array containing the numbers for the outer dofs
static int N_T_RT0_3D_Outer[4] = { 0, 1, 2, 3 };

static char N_T_RT0_3D_String[] = "N_T_RT0_3D";

TFEDesc3D *FE_N_T_RT0_3D_Obj=new TFEDesc3D(N_T_RT0_3D_String, N_T_RT0_3D_NDOF,
                                           N_T_RT0_3D_JointDOF, N_T_RT0_3D_J, 
                                           N_T_RT0_3D_NInner, N_T_RT0_3D_Inner,
                                           N_T_RT0_3D_NOuter, N_T_RT0_3D_Outer);
