// ***********************************************************************
// conforming P4 element with cell bubbles
// ***********************************************************************

// number of degrees of freedom
static int C_T_B4_2D_NDOF = 18;

// number of dofs on the closure of the joints
static int C_T_B4_2D_JointDOF = 5;

// which local dofs are on the joints
static int C_T_B4_2D_J0[5] = { 0, 1,  2,  3, 4 };
static int C_T_B4_2D_J1[5] = { 4, 5,  6,  7, 8 };
static int C_T_B4_2D_J2[5] = { 8, 9, 10, 11, 0 };

static int *C_T_B4_2D_J[3] = { C_T_B4_2D_J0, C_T_B4_2D_J1,
                               C_T_B4_2D_J2 };

// number of inner dofs
static int C_T_B4_2D_NInner = 6;

// array containing the numbers for the inner dofs
static int C_T_B4_2D_Inner[6] = { 12, 13, 14, 15, 16, 17 };

// number of outer dofs
static int C_T_B4_2D_NOuter = 12;

// array containing the numbers for the outer dofs
static int C_T_B4_2D_Outer[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };

static char C_T_B4_2D_String[] = "C_T_B4_2D";

TFEDesc2D *FE_C_T_B4_2D_Obj=new TFEDesc2D(C_T_B4_2D_String, C_T_B4_2D_NDOF,
                              C_T_B4_2D_JointDOF, C_T_B4_2D_J,
                              C_T_B4_2D_NInner, C_T_B4_2D_Inner,
                              C_T_B4_2D_NOuter, C_T_B4_2D_Outer);
