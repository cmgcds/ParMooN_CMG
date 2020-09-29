// ***********************************************************************
// Q3 element, nonconforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int N_Q_Q3_2D_NDOF = 15;

// number of dofs on the closure of the joints
static int N_Q_Q3_2D_JointDOF = 3;

// which local dofs are on the joints
static int N_Q_Q3_2D_J0[3] = { 0, 4,  8 };
static int N_Q_Q3_2D_J1[3] = { 1, 5,  9 };
static int N_Q_Q3_2D_J2[3] = { 2, 6, 10 };
static int N_Q_Q3_2D_J3[3] = { 3, 7, 11 };
 
static int *N_Q_Q3_2D_J[4] = { N_Q_Q3_2D_J0, N_Q_Q3_2D_J1,
                                 N_Q_Q3_2D_J2, N_Q_Q3_2D_J3 };

// number of inner dofs
static int N_Q_Q3_2D_NInner = 3;

// array containing the numbers for the inner dofs (here is no inner dof)
static int N_Q_Q3_2D_Inner[3] = { 12, 13, 14 };

// number of outer dofs
static int N_Q_Q3_2D_NOuter = 12;

// array containing the numbers for the outer dofs
static int N_Q_Q3_2D_Outer[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };

static char N_Q_Q3_2D_String[] = "N_Q_Q3_2D";

TFEDesc2D *FE_N_Q_Q3_2D_Obj=new TFEDesc2D(N_Q_Q3_2D_String, N_Q_Q3_2D_NDOF,
                                        N_Q_Q3_2D_JointDOF, N_Q_Q3_2D_J,
                                        N_Q_Q3_2D_NInner, N_Q_Q3_2D_Inner,
                                        N_Q_Q3_2D_NOuter, N_Q_Q3_2D_Outer);
