//************************************************************
// P5 element with bubbles, conforming, 2D
//************************************************************

// number of degrees of freedom
static int C_T_UL5_2D_NDOF = 30;

// number of dofs on the closure of the joints
static int C_T_UL5_2D_JointDOF = 6;

// which local dofs are on the joints
static int C_T_UL5_2D_J0[6] = {  0,  1,  2,  3,  4,  5 };
static int C_T_UL5_2D_J1[6] = {  5,  6,  7,  8,  9, 10 };
static int C_T_UL5_2D_J2[6] = { 10, 11, 12, 13, 14,  0 };

static int *C_T_UL5_2D_J[3] = { C_T_UL5_2D_J0, C_T_UL5_2D_J1,
                                C_T_UL5_2D_J2 };

// number of inner dofs
static int C_T_UL5_2D_NInner = 15;

// array containing the numbers for the inner dofs
static int C_T_UL5_2D_Inner[15] = { 15, 16, 17, 18, 19, 20, 21, 22, 23,
                                    24, 25, 26, 27, 28, 29 };

// number of outer dofs         
static int C_T_UL5_2D_NOuter = 15;

// array containing the numbers for the inner dofs
static int C_T_UL5_2D_Outer[15] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                                    12, 13, 14 };

static char C_T_UL5_2D_String[] = "C_T_UL5_2D";

TFEDesc2D *FE_C_T_UL5_2D_Obj=new TFEDesc2D(C_T_UL5_2D_String, C_T_UL5_2D_NDOF,
                                C_T_UL5_2D_JointDOF, C_T_UL5_2D_J,
                                C_T_UL5_2D_NInner, C_T_UL5_2D_Inner,
                                C_T_UL5_2D_NOuter, C_T_UL5_2D_Outer);
