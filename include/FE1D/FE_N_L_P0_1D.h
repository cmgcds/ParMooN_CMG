// ***********************************************************************
// P0 element, nonconforming, 1D
// ***********************************************************************

// number of degrees of freedom
static int N_L_P0_1D_NDOF = 1;

// number of dofs on the closure of the joints
static int N_L_P0_1D_JointDOF = 0;

// which local dofs are on the joints
static int *N_L_P0_1D_J[2] = {0, 0};

// number of inner dofs
static int N_L_P0_1D_NInner = 1;

// array containing the numbers for the inner dofs
static int N_L_P0_1D_Inner[1] = { 0 };

static char N_L_P0_1D_String[] = "N_L_P0_1D";

TFEDesc1D *FE_N_L_P0_1D_Obj=new TFEDesc1D(N_L_P0_1D_String, N_L_P0_1D_NDOF, 
                              N_L_P0_1D_JointDOF, N_L_P0_1D_J, 
                              N_L_P0_1D_NInner, N_L_P0_1D_Inner);
