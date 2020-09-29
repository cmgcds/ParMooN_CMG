// ***********************************************************************
// conforming P2 element with cell bubble
// ***********************************************************************

// number of degrees of freedom
static int C_T_B2_2D_NDOF = 7;

// number of dofs on the closure of the joints
static int C_T_B2_2D_JointDOF = 3;

// which local dofs are on the joints
static int C_T_B2_2D_J0[3] = { 0, 1, 2 };
static int C_T_B2_2D_J1[3] = { 2, 4, 5 };
static int C_T_B2_2D_J2[3] = { 5, 3, 0 };

static int *C_T_B2_2D_J[3] = { C_T_B2_2D_J0, C_T_B2_2D_J1,
                              C_T_B2_2D_J2 };

// number of inner dofs
static int C_T_B2_2D_NInner = 1;

// array containing the numbers for the inner dofs
static int C_T_B2_2D_Inner[1] = { 6 };

// number of outer dofs
static int C_T_B2_2D_NOuter = 6;

// array containing the numbers for the outer dofs
static int C_T_B2_2D_Outer[6] = { 0, 1, 2, 3, 4, 5 };

static char C_T_B2_2D_String[] = "C_T_B2_2D";

TFEDesc2D *FE_C_T_B2_2D_Obj=new TFEDesc2D(C_T_B2_2D_String, C_T_B2_2D_NDOF,
                              C_T_B2_2D_JointDOF, C_T_B2_2D_J,
                              C_T_B2_2D_NInner, C_T_B2_2D_Inner,
                              C_T_B2_2D_NOuter, C_T_B2_2D_Outer);
