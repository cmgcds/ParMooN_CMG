// ***********************************************************************
// P1 element, nonconforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int N_T_P1_2D_NDOF = 3;

// number of dofs on the closure of the joints
static int N_T_P1_2D_JointDOF = 1;

// which local dofs are on the joints
static int N_T_P1_2D_J0[1] = { 0 };
static int N_T_P1_2D_J1[1] = { 1 };
static int N_T_P1_2D_J2[1] = { 2 };

static int *N_T_P1_2D_J[3] = { N_T_P1_2D_J0, N_T_P1_2D_J1,
                                 N_T_P1_2D_J2 };

// number of inner dofs
static int N_T_P1_2D_NInner = 0;

// array containing the numbers for the inner dofs (here is no inner dof)
static int *N_T_P1_2D_Inner = NULL;

// number of outer dofs
static int N_T_P1_2D_NOuter = 3;

// array containing the numbers for the outer dofs
static int N_T_P1_2D_Outer[3] = { 0, 1, 2 };

static char N_T_P1_2D_String[] = "N_T_P1_2D";

TFEDesc2D *FE_N_T_P1_2D_Obj=new TFEDesc2D(N_T_P1_2D_String, N_T_P1_2D_NDOF,
                                        N_T_P1_2D_JointDOF, N_T_P1_2D_J,
                                        N_T_P1_2D_NInner, N_T_P1_2D_Inner,
                                        N_T_P1_2D_NOuter, N_T_P1_2D_Outer);
