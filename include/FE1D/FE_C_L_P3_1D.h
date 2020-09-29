// ***********************************************************************
// P3 element, conforming, 1D
// ***********************************************************************

// number of degrees of freedom
static int C_L_P3_1D_NDOF = 4;

// number of dofs on the closure of the joints
static int C_L_P3_1D_JointDOF = 1;

// which local dofs are on the joints
static int C_L_P3_1D_J0[1] = {  0 };
static int C_L_P3_1D_J1[1] = {  3 };

static int *C_L_P3_1D_J[2] = { C_L_P3_1D_J0, C_L_P3_1D_J1 };

// number of inner dofs
static int C_L_P3_1D_NInner = 2;

// array containing the numbers for the inner dofs
static int C_L_P3_1D_Inner[2] = { 1, 2 };

static char C_L_P3_1D_String[] = "C_L_P3_1D";

TFEDesc1D *FE_C_L_P3_1D_Obj=new TFEDesc1D(C_L_P3_1D_String, C_L_P3_1D_NDOF, 
                              C_L_P3_1D_JointDOF, C_L_P3_1D_J, 
                              C_L_P3_1D_NInner, C_L_P3_1D_Inner);
