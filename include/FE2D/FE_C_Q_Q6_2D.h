//************************************************************
// Q6 element, conforming, 2D
//************************************************************

// number of degrees of freedom
static int C_Q_Q6_2D_NDOF = 49;

// number of dofs on the closure of the joints
static int C_Q_Q6_2D_JointDOF = 7;

// which local dofs are on the joints
static int C_Q_Q6_2D_J0[7] = { 0, 1, 2, 3, 4, 5, 6 };
static int C_Q_Q6_2D_J1[7] = { 6, 13, 20, 27, 34, 41, 48 };
static int C_Q_Q6_2D_J2[7] = { 48, 47, 46, 45, 44, 43, 42 };
static int C_Q_Q6_2D_J3[7] = { 42, 35, 28, 21, 14, 7, 0 };

static int *C_Q_Q6_2D_J[7] = {  C_Q_Q6_2D_J0,  C_Q_Q6_2D_J1,
                              C_Q_Q6_2D_J2,  C_Q_Q6_2D_J3 };

// number of inner dofs
static int C_Q_Q6_2D_NInner = 25;

// array containing the numbers for the inner dofs
static int C_Q_Q6_2D_Inner[25] = { 8, 9, 10, 11, 12, 15, 16, 17, 18, 19, 22,
                                   23, 24, 25, 26, 29, 30, 31, 32, 33, 36,
                                   37, 38, 39, 40 };

// number of outer dofs
static int C_Q_Q6_2D_NOuter = 24;

// array containing the numbers for the outer dofs
static int C_Q_Q6_2D_Outer[24] = { 0, 1, 2, 3, 4, 5, 6, 7, 13, 14, 20, 21, 
                                   27, 28, 34, 35, 41, 42, 43, 44, 45, 46,
                                   47, 48 };

static char C_Q_Q6_2D_String[] = "C_Q_Q6_2D";

TFEDesc2D *FE_C_Q_Q6_2D_Obj=new TFEDesc2D(C_Q_Q6_2D_String, C_Q_Q6_2D_NDOF,
                                        C_Q_Q6_2D_JointDOF, C_Q_Q6_2D_J,
                                        C_Q_Q6_2D_NInner, C_Q_Q6_2D_Inner,
                                        C_Q_Q6_2D_NOuter, C_Q_Q6_2D_Outer);
