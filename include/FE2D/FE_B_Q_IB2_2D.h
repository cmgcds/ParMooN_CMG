// ***********************************************************************
// internal bubble of degree 2 (in the sense of Q2)
// ***********************************************************************

// number of degrees of freedom
static int B_Q_IB2_2D_NDOF = 1;

// number of dofs on the closure of the joints
static int B_Q_IB2_2D_JointDOF = 0;

// which local dofs are on the joints
static int *B_Q_IB2_2D_J[3] = { NULL, NULL, NULL };

// number of inner dofs
static int B_Q_IB2_2D_NInner = 1;

// array containing the numbers for the inner dofs 
static int B_Q_IB2_2D_Inner[1] = { 0 };

static char B_Q_IB2_2D_String[] = "B_Q_IB2_2D";

TFEDesc2D *FE_B_Q_IB2_2D_Obj=new TFEDesc2D(B_Q_IB2_2D_String, B_Q_IB2_2D_NDOF,
                                           B_Q_IB2_2D_JointDOF, B_Q_IB2_2D_J,
                                           B_Q_IB2_2D_NInner, B_Q_IB2_2D_Inner);
