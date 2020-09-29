// ***********************************************************************
// P2 element, conforming, 3D
// ***********************************************************************

// number of degrees of freedom
static int C_T_P2_3D_NDOF = 10;

// number of dofs on the closure of the joints
static int C_T_P2_3D_JointDOF = 6;

// which local dofs are on the joints
static int C_T_P2_3D_J0[6] = {  0,  1,  2,  3,  4,  5 };
static int C_T_P2_3D_J1[6] = {  0,  6,  9,  1,  7,  2 };
static int C_T_P2_3D_J2[6] = {  5,  4,  2,  8,  7,  9 };
static int C_T_P2_3D_J3[6] = {  0,  3,  5,  6,  8,  9 };

static int *C_T_P2_3D_J[4] = { C_T_P2_3D_J0, C_T_P2_3D_J1,
                             C_T_P2_3D_J2, C_T_P2_3D_J3 };


#ifdef _MPI   
// number of dofs on the closure of the edges
static int C_T_P2_3D_EdgeDOF = 3;

// which local dofs are on the joints
static int C_T_P2_3D_E0[3] = { 0, 1, 2 };
static int C_T_P2_3D_E1[3] = { 2, 4, 5 };
static int C_T_P2_3D_E2[3] = { 5, 3, 0 };
static int C_T_P2_3D_E3[3] = { 0, 6, 9 };
static int C_T_P2_3D_E4[3] = { 2, 7, 9 };
static int C_T_P2_3D_E5[3] = { 5, 8, 9 };


static int *C_T_P2_3D_E[6] = { C_T_P2_3D_E0, C_T_P2_3D_E1, C_T_P2_3D_E2, C_T_P2_3D_E3,
                               C_T_P2_3D_E4, C_T_P2_3D_E5};

// number of dofs on the closure of the vertices
static int C_T_P2_3D_VertDOF = 1;

// array containing the numbers for the vertices dofs
static int C_T_P2_3D_Vert[4] =  {0, 2, 5, 9};
#endif



// number of inner dofs
static int C_T_P2_3D_NInner = 0;

// array containing the numbers for the inner dofs (here is no inner dof)
static int *C_T_P2_3D_Inner = NULL;

static char C_T_P2_3D_String[] = "C_T_P2_3D";

TFEDesc3D *FE_C_T_P2_3D_Obj=new TFEDesc3D(C_T_P2_3D_String, C_T_P2_3D_NDOF, 
                                C_T_P2_3D_JointDOF,
                                C_T_P2_3D_J, C_T_P2_3D_NInner, C_T_P2_3D_Inner
#ifdef _MPI
                                ,C_T_P2_3D_EdgeDOF,  C_T_P2_3D_E, C_T_P2_3D_VertDOF,
                                 C_T_P2_3D_Vert
#endif
                                 );


