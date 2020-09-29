//************************************************************
// Q5 element, conforming, 2D
//************************************************************

// number of degrees of freedom
static int C_Q_Q5_2D_NDOF = 36;

// number of dofs on the closure of the joints
static int C_Q_Q5_2D_JointDOF = 6;

// which local dofs are on the joints
static int C_Q_Q5_2D_J0[6] = { 0, 1, 2, 3, 4, 5 };
static int C_Q_Q5_2D_J1[6] = { 5, 11, 17, 23, 29, 35 };
static int C_Q_Q5_2D_J2[6] = { 35, 34, 33, 32, 31, 30 };
static int C_Q_Q5_2D_J3[6] = { 30, 24, 18, 12, 6, 0 };

static int *C_Q_Q5_2D_J[6] = {  C_Q_Q5_2D_J0,  C_Q_Q5_2D_J1,
                              C_Q_Q5_2D_J2,  C_Q_Q5_2D_J3 };

// number of inner dofs
static int C_Q_Q5_2D_NInner = 16;

// array containing the numbers for the inner dofs
static int C_Q_Q5_2D_Inner[16] = { 7, 8, 9, 10, 13, 14, 15, 16, 19, 20,
                                   21, 22, 25, 26, 27, 28 };

// number of outer dofs         
static int C_Q_Q5_2D_NOuter = 20;

// array containing the numbers for the inner dofs
static int C_Q_Q5_2D_Outer[20] = { 0, 1, 2, 3, 4, 5, 6, 11, 12, 17, 18, 
                                   23, 24, 29, 30, 31, 32, 33, 34, 35 }; 
                                  

static char C_Q_Q5_2D_String[] = "C_Q_Q5_2D";

TFEDesc2D *FE_C_Q_Q5_2D_Obj=new TFEDesc2D(C_Q_Q5_2D_String, C_Q_Q5_2D_NDOF,
                                C_Q_Q5_2D_JointDOF, C_Q_Q5_2D_J,
                                C_Q_Q5_2D_NInner, C_Q_Q5_2D_Inner,
                                C_Q_Q5_2D_NOuter, C_Q_Q5_2D_Outer);
