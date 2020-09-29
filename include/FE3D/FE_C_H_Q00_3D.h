// ***********************************************************************
// Q00 element, discontinuous, 3D
// ***********************************************************************

// number of degrees of freedom
static int C_H_Q00_3D_NDOF = 1;

// number of dofs on the closure of the joints
static int C_H_Q00_3D_JointDOF = 0;

// which local dofs are on the joints
static int *C_H_Q00_3D_J0 = NULL;
static int *C_H_Q00_3D_J1 = NULL;
static int *C_H_Q00_3D_J2 = NULL;
static int *C_H_Q00_3D_J3 = NULL;
static int *C_H_Q00_3D_J4 = NULL;
static int *C_H_Q00_3D_J5 = NULL;

static int *C_H_Q00_3D_J[6] = { C_H_Q00_3D_J0, C_H_Q00_3D_J1,
                               C_H_Q00_3D_J2, C_H_Q00_3D_J3,
                               C_H_Q00_3D_J4, C_H_Q00_3D_J5};

// number of inner dofs
static int C_H_Q00_3D_NInner = 1;

// array containing the numbers for the inner dofs
static int C_H_Q00_3D_Inner[1] = { 0 };

static char C_H_Q00_3D_String[] = "C_H_Q00_3D";

TFEDesc3D *FE_C_H_Q00_3D_Obj=new TFEDesc3D(C_H_Q00_3D_String, C_H_Q00_3D_NDOF, 
                                C_H_Q00_3D_JointDOF,
                                C_H_Q00_3D_J, C_H_Q00_3D_NInner, C_H_Q00_3D_Inner);
