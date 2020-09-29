// ***********************************************************************
// Brezzi-Douglas-Duran-Fortin element of first order on hexahedra, 3D
// ***********************************************************************

// number of degrees of freedom
static int N_H_BDDF1_3D_NDOF = 18;

// number of dofs on the closure of each joints
static int N_H_BDDF1_3D_JointDOF = 3;

// which local dofs are on the joints
static int N_H_BDDF1_3D_J0[4] = { 0, 1, 2 };
static int N_H_BDDF1_3D_J1[4] = { 3, 4, 5 };
static int N_H_BDDF1_3D_J2[4] = { 6, 7, 8 };
static int N_H_BDDF1_3D_J3[4] = { 9,10,11 };
static int N_H_BDDF1_3D_J4[4] = {12,13,14 };
static int N_H_BDDF1_3D_J5[4] = {15,16,17 };

static int *N_H_BDDF1_3D_J[6] = { N_H_BDDF1_3D_J0, N_H_BDDF1_3D_J1,
                                N_H_BDDF1_3D_J2, N_H_BDDF1_3D_J3,
                                N_H_BDDF1_3D_J4, N_H_BDDF1_3D_J5 };

// number of inner dofs
static int N_H_BDDF1_3D_NInner = 0;

// array containing the numbers for the inner dofs
static int *N_H_BDDF1_3D_Inner = NULL;

// number of outer dofs
static int N_H_BDDF1_3D_NOuter = 18;

// array containing the numbers for the outer dofs
static int N_H_BDDF1_3D_Outer[24] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17 };

static char N_H_BDDF1_3D_String[] = "N_H_BDDF1_3D";

TFEDesc3D *FE_N_H_BDDF1_3D_Obj=new TFEDesc3D(N_H_BDDF1_3D_String, N_H_BDDF1_3D_NDOF,
                                           N_H_BDDF1_3D_JointDOF, N_H_BDDF1_3D_J,
                                           N_H_BDDF1_3D_NInner, N_H_BDDF1_3D_Inner,
                                           N_H_BDDF1_3D_NOuter, N_H_BDDF1_3D_Outer);
