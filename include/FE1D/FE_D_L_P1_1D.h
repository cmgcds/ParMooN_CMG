// ***********************************************************************
// P1 discontinuous element, 1D
// ***********************************************************************

// number of degrees of freedom
static int D_L_P1_1D_NDOF = 2;

// number of dofs on the closure of the joints
static int D_L_P1_1D_JointDOF = 0;

// which local dofs are on the joints
static int *D_L_P1_1D_J0 = NULL;
static int *D_L_P1_1D_J1 = NULL;

static int *D_L_P1_1D_J[2] = { D_L_P1_1D_J0, D_L_P1_1D_J1 };

// number of inner dofs
static int D_L_P1_1D_NInner = 2;

// array containing the numbers for the inner dofs
static int D_L_P1_1D_Inner[2] = {0, 1};

static char D_L_P1_1D_String[] = "D_L_P1_1D";

TFEDesc1D *FE_D_L_P1_1D_Obj=new TFEDesc1D(D_L_P1_1D_String, D_L_P1_1D_NDOF, 
                              D_L_P1_1D_JointDOF, D_L_P1_1D_J, 
                              D_L_P1_1D_NInner, D_L_P1_1D_Inner);
