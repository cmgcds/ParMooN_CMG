// ***********************************************************************
// conforming P6 element
// ***********************************************************************

// number of degrees of freedom
static int C_T_P6_2D_NDOF = 28;

// number of dofs on the closure of the joints
static int C_T_P6_2D_JointDOF = 7;

// which local dofs are on the joints
static int C_T_P6_2D_J0[7] = {  0,  1,  2,  3,  4,  5,  6 };
static int C_T_P6_2D_J1[7] = {  6, 12, 17, 21, 24, 26, 27 };
static int C_T_P6_2D_J2[7] = { 27, 25, 22, 18, 13,  7,  0 };

static int *C_T_P6_2D_J[3] = { C_T_P6_2D_J0, C_T_P6_2D_J1,
                               C_T_P6_2D_J2 };

// number of inner dofs
static int C_T_P6_2D_NInner = 10;

// array containing the numbers for the inner dofs
static int C_T_P6_2D_Inner[10] = { 8, 9, 10, 11, 14, 15, 16, 19, 20, 23 };

// number of outer dofs
static int C_T_P6_2D_NOuter = 18;

// array containing the numbers for the outer dofs
static int C_T_P6_2D_Outer[18] = { 0,  1,  2,  3,  4,  5,
                                   6, 12, 17, 21, 24, 26,
                                  27, 25, 18, 13,  7,  0 }; 

static char C_T_P6_2D_String[] = "C_T_P6_2D";

TFEDesc2D *FE_C_T_P6_2D_Obj=new TFEDesc2D(C_T_P6_2D_String, C_T_P6_2D_NDOF,
                              C_T_P6_2D_JointDOF, C_T_P6_2D_J,
                              C_T_P6_2D_NInner, C_T_P6_2D_Inner,
                              C_T_P6_2D_NOuter, C_T_P6_2D_Outer);
