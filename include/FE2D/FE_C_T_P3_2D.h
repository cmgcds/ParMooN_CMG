// ***********************************************************************
// P3 element, conforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int C_T_P3_2D_NDOF = 10;

// number of dofs on the closure of the joints
static int C_T_P3_2D_JointDOF = 4;

// which local dofs are on the joints
static int C_T_P3_2D_J0[4] = { 0, 1, 2, 3 };
static int C_T_P3_2D_J1[4] = { 3, 6, 8, 9 };
static int C_T_P3_2D_J2[4] = { 9, 7, 4, 0 };

static int *C_T_P3_2D_J[3] = { C_T_P3_2D_J0, C_T_P3_2D_J1,
                              C_T_P3_2D_J2 };

// number of inner dofs
static int C_T_P3_2D_NInner = 1;

// array containing the numbers for the inner dofs
static int C_T_P3_2D_Inner[1] = { 5 };

static char C_T_P3_2D_String[] = "C_T_P3_2D";

TFEDesc2D *FE_C_T_P3_2D_Obj=new TFEDesc2D(C_T_P3_2D_String, C_T_P3_2D_NDOF, C_T_P3_2D_JointDOF,
                              C_T_P3_2D_J, C_T_P3_2D_NInner, C_T_P3_2D_Inner);
