// ***********************************************************************
// Q1 element, nonconforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int N_Q_Q1_2D_NDOF = 4;

// number of dofs on the closure of the joints
static int N_Q_Q1_2D_JointDOF = 1;

// which local dofs are on the joints
static int N_Q_Q1_2D_J0[1] = { 0 };
static int N_Q_Q1_2D_J1[1] = { 1 };
static int N_Q_Q1_2D_J2[1] = { 2 };
static int N_Q_Q1_2D_J3[1] = { 3 };
 
static int *N_Q_Q1_2D_J[4] = { N_Q_Q1_2D_J0, N_Q_Q1_2D_J1,
                                 N_Q_Q1_2D_J2, N_Q_Q1_2D_J3 };

// number of inner dofs
static int N_Q_Q1_2D_NInner = 0;

// array containing the numbers for the inner dofs (here is no inner dof)
static int *N_Q_Q1_2D_Inner = NULL;

// number of outer dofs
static int N_Q_Q1_2D_NOuter = 4;

// array containing the numbers for the outer dofs
static int N_Q_Q1_2D_Outer[4] = { 0, 1, 2, 3 };

static char N_Q_Q1_2D_String[] = "N_Q_Q1_2D";

TFEDesc2D *FE_N_Q_Q1_2D_Obj=new TFEDesc2D(N_Q_Q1_2D_String, N_Q_Q1_2D_NDOF,
                                        N_Q_Q1_2D_JointDOF, N_Q_Q1_2D_J,
                                        N_Q_Q1_2D_NInner, N_Q_Q1_2D_Inner,
                                        N_Q_Q1_2D_NOuter, N_Q_Q1_2D_Outer);
