// ***********************************************************************
// Brezzi-Douglas-Duran-Fortin element of second order on tetrahedra, 3D
// ***********************************************************************

// number of degrees of freedom
static int N_T_BDDF2_3D_NDOF = 30;

// number of dofs on the closure of each joints
static int N_T_BDDF2_3D_JointDOF = 6;

// which local dofs are on the joints
static int N_T_BDDF2_3D_J0[6] = { 0, 1, 2, 3, 4, 5};
static int N_T_BDDF2_3D_J1[6] = { 6, 7, 8, 9,10,11};
static int N_T_BDDF2_3D_J2[6] = {12,13,14,15,16,17};
static int N_T_BDDF2_3D_J3[6] = {18,19,20,21,22,23};

static int *N_T_BDDF2_3D_J[4] = { N_T_BDDF2_3D_J0, N_T_BDDF2_3D_J1,
                                N_T_BDDF2_3D_J2, N_T_BDDF2_3D_J3 };

// number of inner dofs
static int N_T_BDDF2_3D_NInner = 6;

// array containing the numbers for the inner dofs (here is no inner dof)
static int N_T_BDDF2_3D_Inner[6] = {24,25,26,27,28,29};

// number of outer dofs
static int N_T_BDDF2_3D_NOuter = 24;

// array containing the numbers for the outer dofs
static int N_T_BDDF2_3D_Outer[24] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23 };

static char N_T_BDDF2_3D_String[] = "N_T_BDDF2_3D";

TFEDesc3D *FE_N_T_BDDF2_3D_Obj=new TFEDesc3D(N_T_BDDF2_3D_String, N_T_BDDF2_3D_NDOF,
                                           N_T_BDDF2_3D_JointDOF, N_T_BDDF2_3D_J,
                                           N_T_BDDF2_3D_NInner, N_T_BDDF2_3D_Inner,
                                           N_T_BDDF2_3D_NOuter, N_T_BDDF2_3D_Outer);
