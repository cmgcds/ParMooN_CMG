// ***********************************************************************
// P2 element, conforming, 2D
// ***********************************************************************

// number of degrees of freedom
static int C_T_P2_2D_NDOF = 6;

// number fo dofs on the closure of joints
static int C_T_P2_2D_JointDOF = 3;

// which local dofs are on the joints
static int C_T_P2_2D_J0[3]={ 0, 1, 2 };
static int C_T_P2_2D_J1[3]={ 2, 4, 5 };
static int C_T_P2_2D_J2[3]={ 5, 3, 0 };

static int *C_T_P2_2D_J[3]={ C_T_P2_2D_J0, C_T_P2_2D_J1,
                           C_T_P2_2D_J2 };

// number of inner dofs
static int C_T_P2_2D_NInner = 0;

// array containing the numbers for the inner dofs
static int *C_T_P2_2D_Inner = NULL;

// array containing the numbers for the inner dofs
static char C_T_P2_2D_String[] = "C_T_P2_2D";

TFEDesc2D *FE_C_T_P2_2D_Obj=new TFEDesc2D(C_T_P2_2D_String, C_T_P2_2D_NDOF, C_T_P2_2D_JointDOF,
                              C_T_P2_2D_J, C_T_P2_2D_NInner, C_T_P2_2D_Inner);

