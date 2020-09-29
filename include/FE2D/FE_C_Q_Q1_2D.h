// ***********************************************************************
// Q1 element, conforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int C_Q_Q1_2D_NDOF = 4;

// number of dofs on the closure of the joints
static int C_Q_Q1_2D_JointDOF = 2;

// which local dofs are on the joints
static int C_Q_Q1_2D_J0[2] = { 0, 1 };
static int C_Q_Q1_2D_J1[2] = { 1, 3 };
static int C_Q_Q1_2D_J2[2] = { 3, 2 };
static int C_Q_Q1_2D_J3[2] = { 2, 0 };

static int *C_Q_Q1_2D_J[4] = { C_Q_Q1_2D_J0, C_Q_Q1_2D_J1,
                              C_Q_Q1_2D_J2, C_Q_Q1_2D_J3 };

// number of inner dofs
static int C_Q_Q1_2D_NInner = 0;

// array containing the numbers for the inner dofs (here is no inner dof)
static int *C_Q_Q1_2D_Inner = NULL;

// number of outer dofs
static int C_Q_Q1_2D_NOuter = 4;

// array containing the numbers for the outer dofs
static int C_Q_Q1_2D_Outer[4] = { 0, 1, 2, 3 };

static char C_Q_Q1_2D_String[] = "C_Q_Q1_2D";

TFEDesc2D *FE_C_Q_Q1_2D_Obj=new TFEDesc2D(C_Q_Q1_2D_String, C_Q_Q1_2D_NDOF,
                                C_Q_Q1_2D_JointDOF, C_Q_Q1_2D_J,
                                C_Q_Q1_2D_NInner, C_Q_Q1_2D_Inner,
                                C_Q_Q1_2D_NOuter, C_Q_Q1_2D_Outer);
