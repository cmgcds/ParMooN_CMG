// ***********************************************************************
// UL2S element, conforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int C_Q_M2_2D_NDOF = 8;

// number of dofs on the closure of the joints
static int C_Q_M2_2D_JointDOF = 3;

// which local dofs are on the joints
static int C_Q_M2_2D_J0[3] = { 0, 1, 2 };
static int C_Q_M2_2D_J1[3] = { 2, 3, 4 };
static int C_Q_M2_2D_J2[3] = { 4, 5, 6 };
static int C_Q_M2_2D_J3[3] = { 6, 7, 0 };

static int *C_Q_M2_2D_J[4] = { C_Q_M2_2D_J0, C_Q_M2_2D_J1,
                               C_Q_M2_2D_J2, C_Q_M2_2D_J3 };

// number of inner dofs
static int C_Q_M2_2D_NInner = 0;

// array containing the numbers for the inner dofs
static int *C_Q_M2_2D_Inner = NULL;

// number of outer dofs
static int C_Q_M2_2D_NOuter = 8;

// array containing the numbers for the outer dofs
static int C_Q_M2_2D_Outer[8] = { 0, 1, 2, 3, 5, 6, 7, 8 };

static char C_Q_M2_2D_String[] = "C_Q_M2_2D";

TFEDesc2D *FE_C_Q_M2_2D_Obj=new TFEDesc2D(C_Q_M2_2D_String, C_Q_M2_2D_NDOF,
                                C_Q_M2_2D_JointDOF, C_Q_M2_2D_J,
                                C_Q_M2_2D_NInner, C_Q_M2_2D_Inner,
                                C_Q_M2_2D_NOuter, C_Q_M2_2D_Outer);
