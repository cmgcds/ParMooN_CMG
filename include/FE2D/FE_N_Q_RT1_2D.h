// ***********************************************************************
// Q1 Raviart-Thomas vector element, nonconforming , 2D
// History:  17.11.2011 implementation (Ulrich)
// ***********************************************************************

// number of degrees of freedom
static int N_Q_RT1_2D_NDOF = 12;

// number of dofs on the closure of the joints
static int N_Q_RT1_2D_JointDOF = 2;

// which local dofs are on the joints
static int N_Q_RT1_2D_J0[2] = { 0, 1 };
static int N_Q_RT1_2D_J1[2] = { 2, 3 };
static int N_Q_RT1_2D_J2[2] = { 4, 5 };
static int N_Q_RT1_2D_J3[2] = { 6, 7 };
 
static int *N_Q_RT1_2D_J[4] = { N_Q_RT1_2D_J0, N_Q_RT1_2D_J1,
                                 N_Q_RT1_2D_J2, N_Q_RT1_2D_J3 };

// number of inner dofs
static int N_Q_RT1_2D_NInner = 4;

// array containing the numbers for the inner dofs 
static int N_Q_RT1_2D_Inner[4] = {8,9,10,11};

// number of outer dofs
static int N_Q_RT1_2D_NOuter = 8;

// array containing the numbers for the outer dofs
static int N_Q_RT1_2D_Outer[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };

static char N_Q_RT1_2D_String[] = "N_Q_RT1_2D";

TFEDesc2D *FE_N_Q_RT1_2D_Obj=new TFEDesc2D(N_Q_RT1_2D_String, N_Q_RT1_2D_NDOF,
                                        N_Q_RT1_2D_JointDOF, N_Q_RT1_2D_J,
                                        N_Q_RT1_2D_NInner, N_Q_RT1_2D_Inner,
                                        N_Q_RT1_2D_NOuter, N_Q_RT1_2D_Outer);
