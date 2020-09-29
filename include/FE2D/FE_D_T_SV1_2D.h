// ***********************************************************************
// P1 element, discontinuous, 2D, on a subdivision into three triangles
// ***********************************************************************

// number of degrees of freedom
static int D_T_SV1_2D_NDOF = 9;

// number fo dofs on the closure of joints
static int D_T_SV1_2D_JointDOF = 0;

// which local dofs are on the joints
static int *D_T_SV1_2D_J0 = NULL;
static int *D_T_SV1_2D_J1 = NULL;
static int *D_T_SV1_2D_J2 = NULL;

static int *D_T_SV1_2D_J[3]={ D_T_SV1_2D_J0, D_T_SV1_2D_J1,
                              D_T_SV1_2D_J2 };

// number of inner dofs
static int D_T_SV1_2D_NInner = 9;

// array containing the numbers for the inner dofs
static int D_T_SV1_2D_Inner[9] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };

// describing string
static char D_T_SV1_2D_String[] = "D_T_SV1_2D";

TFEDesc2D *FE_D_T_SV1_2D_Obj=new TFEDesc2D(
                D_T_SV1_2D_String, D_T_SV1_2D_NDOF, D_T_SV1_2D_JointDOF,
                D_T_SV1_2D_J, D_T_SV1_2D_NInner, D_T_SV1_2D_Inner);

