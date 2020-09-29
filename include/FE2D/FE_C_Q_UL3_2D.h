// ***********************************************************************
// Q3 element with bubbles, conforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int C_Q_UL3_2D_NDOF = 18;

// number of dofs on the closure of the joints
static int C_Q_UL3_2D_JointDOF = 4;

// which local dofs are on the joints
static int C_Q_UL3_2D_J0[4] = {  0,  1,  2,  3 };
static int C_Q_UL3_2D_J1[4] = {  3,  4,  5,  6 };
static int C_Q_UL3_2D_J2[4] = {  6,  7,  8,  9 };
static int C_Q_UL3_2D_J3[4] = {  9, 10, 11,  0 };

static int *C_Q_UL3_2D_J[4] = { C_Q_UL3_2D_J0, C_Q_UL3_2D_J1, 
                                C_Q_UL3_2D_J2, C_Q_UL3_2D_J3 };

// number of inner dofs
static int C_Q_UL3_2D_NInner = 6;

// array containing the numbers for the inner dofs
static int C_Q_UL3_2D_Inner[6] = { 12, 13, 14, 15, 16, 17 };

// number of outer dofs
static int C_Q_UL3_2D_NOuter = 12;

// array containing the numbers for the outer dofs
static int C_Q_UL3_2D_Outer[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };

static char C_Q_UL3_2D_String[] = "C_Q_UL3_2D";

TFEDesc2D *FE_C_Q_UL3_2D_Obj=new TFEDesc2D(C_Q_UL3_2D_String, C_Q_UL3_2D_NDOF,
                                C_Q_UL3_2D_JointDOF, C_Q_UL3_2D_J,
                                C_Q_UL3_2D_NInner, C_Q_UL3_2D_Inner,
                                C_Q_UL3_2D_NOuter, C_Q_UL3_2D_Outer);
