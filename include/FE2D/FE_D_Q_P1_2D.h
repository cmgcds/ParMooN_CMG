// ***********************************************************************
// P1 element, discontinous, 2D, quadrilateral
// ***********************************************************************

// number of degrees of freedom
static int D_Q_P1_2D_NDOF = 3;

// number of dofs on the closure of the joints
static int D_Q_P1_2D_JointDOF = 0;

// which local dofs are on the joints
static int *D_Q_P1_2D_J0 = NULL;
static int *D_Q_P1_2D_J1 = NULL;
static int *D_Q_P1_2D_J2 = NULL;
static int *D_Q_P1_2D_J3 = NULL;

static int *D_Q_P1_2D_J[4] = { D_Q_P1_2D_J0, D_Q_P1_2D_J1,
                             D_Q_P1_2D_J2, D_Q_P1_2D_J3 };

// number of inner dofs
static int D_Q_P1_2D_NInner = 3;

// array containing the numbers for the inner dofs (here is no inner dof)
static int D_Q_P1_2D_Inner[3] = { 0, 1, 2 };

static char D_Q_P1_2D_String[] = "D_Q_P1_2D";

TFEDesc2D *FE_D_Q_P1_2D_Obj=new TFEDesc2D(D_Q_P1_2D_String, D_Q_P1_2D_NDOF, D_Q_P1_2D_JointDOF,
                              D_Q_P1_2D_J, D_Q_P1_2D_NInner, D_Q_P1_2D_Inner);
