// ***********************************************************************
// Q2 Raviart-Thomas vector element, nonconforming , 2D
// History:  05.12.2011 implementation (Ulrich)
// ***********************************************************************

// number of degrees of freedom
static int N_Q_RT3_2D_NDOF = 40;

// number of dofs on the closure of the joints
static int N_Q_RT3_2D_JointDOF = 4;

// which local dofs are on the joints
static int N_Q_RT3_2D_J0[4] = { 0, 1, 2, 3 };
static int N_Q_RT3_2D_J1[4] = { 4, 5, 6, 7 };
static int N_Q_RT3_2D_J2[4] = { 8, 9, 10,11};
static int N_Q_RT3_2D_J3[4] = { 12,13,14,15};

 
static int *N_Q_RT3_2D_J[4] = { N_Q_RT3_2D_J0,
                                N_Q_RT3_2D_J1,
                                N_Q_RT3_2D_J2,
                                N_Q_RT3_2D_J3
                              };

// number of inner dofs
static int N_Q_RT3_2D_NInner = 24;

// array containing the numbers for the inner dofs
static int N_Q_RT3_2D_Inner[24] = {16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39};

// number of outer dofs (dofs on edges)
static int N_Q_RT3_2D_NOuter = 16;

// array containing the numbers for the outer dofs
static int N_Q_RT3_2D_Outer[16] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};

static char N_Q_RT3_2D_String[] = "N_Q_RT3_2D";

TFEDesc2D *FE_N_Q_RT3_2D_Obj=new TFEDesc2D(N_Q_RT3_2D_String, N_Q_RT3_2D_NDOF,
                                        N_Q_RT3_2D_JointDOF, N_Q_RT3_2D_J,
                                        N_Q_RT3_2D_NInner, N_Q_RT3_2D_Inner,
                                        N_Q_RT3_2D_NOuter, N_Q_RT3_2D_Outer);
