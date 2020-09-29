// ***********************************************************************
// Q4 element, conforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int C_Q_Q4_2D_NDOF = 25;

// number of dofs on the closure of the joints
static int C_Q_Q4_2D_JointDOF = 5;

// which local dofs are on the joints
static int C_Q_Q4_2D_J0[5] = {  0,  1,  2,  3,  4 };
static int C_Q_Q4_2D_J1[5] = {  4,  9, 14, 19, 24 };
static int C_Q_Q4_2D_J2[5] = { 24, 23, 22, 21, 20 };
static int C_Q_Q4_2D_J3[5] = { 20, 15, 10,  5,  0 };

static int *C_Q_Q4_2D_J[4] = { C_Q_Q4_2D_J0, C_Q_Q4_2D_J1, 
                             C_Q_Q4_2D_J2, C_Q_Q4_2D_J3 };

// number of inner dofs
static int C_Q_Q4_2D_NInner = 9;

// array containing the numbers for the inner dofs
static int C_Q_Q4_2D_Inner[9] = { 6, 7, 8, 11, 12, 13, 16, 17, 18 };

// number of outer dofs
static int C_Q_Q4_2D_NOuter = 16;

// array containing the numbers for the outer dofs
static int C_Q_Q4_2D_Outer[16] = { 0, 1, 2, 3, 4, 5, 9, 10, 14, 15, 19,
                                   20, 21, 22, 23, 24 };

static char C_Q_Q4_2D_String[] = "C_Q_Q4_2D";

TFEDesc2D *FE_C_Q_Q4_2D_Obj=new TFEDesc2D(C_Q_Q4_2D_String, C_Q_Q4_2D_NDOF,
                                C_Q_Q4_2D_JointDOF, C_Q_Q4_2D_J,
                                C_Q_Q4_2D_NInner, C_Q_Q4_2D_Inner,
                                C_Q_Q4_2D_NOuter, C_Q_Q4_2D_Outer);
