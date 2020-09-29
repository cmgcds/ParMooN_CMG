// ***********************************************************************
// Q3 element, conforming, 3D
// ***********************************************************************

// number of degrees of freedom
static int C_H_Q3_3D_NDOF = 64;

// number of dofs on the closure of the joints
static int C_H_Q3_3D_JointDOF = 16;

// which local dofs are on the joints
static int C_H_Q3_3D_J0[16] = {  0,  1,  2,  3,  4,  5,  6,  7,
                               8,  9, 10, 11, 12, 13, 14, 15 };
static int C_H_Q3_3D_J1[16] = {  0, 16, 32, 48,  1, 17, 33, 49,
                               2, 18, 34, 50,  3, 19, 35, 51 };
static int C_H_Q3_3D_J2[16] = {  3, 19, 35, 51,  7, 23, 39, 55,
                              11, 27, 43, 59, 15, 31, 47, 63 };
static int C_H_Q3_3D_J3[16] = { 15, 31, 47, 63, 14, 30, 46, 62,
                              13, 29, 45, 61, 12, 28, 44, 60 };
static int C_H_Q3_3D_J4[16] = {  0,  4,  8, 12, 16, 20, 24, 28,
                              32, 36, 40, 44, 48, 52, 56, 60 };
static int C_H_Q3_3D_J5[16] = { 48, 52, 56, 60, 49, 53, 57, 61,
                              50, 54, 58, 62, 51, 55, 59, 63 };

static int *C_H_Q3_3D_J[6] = { C_H_Q3_3D_J0, C_H_Q3_3D_J1,
                             C_H_Q3_3D_J2, C_H_Q3_3D_J3,
                             C_H_Q3_3D_J4, C_H_Q3_3D_J5};

// number of inner dofs
static int C_H_Q3_3D_NInner = 8;

// array containing the numbers for the inner dofs (here is no inner dof)
static int C_H_Q3_3D_Inner[8] = { 21, 22, 25, 26, 37, 38, 41, 42 };

// number of outer dof
static int C_H_Q3_3D_NOuter = 56;

// array containing the numbers for the outer dofs
static int C_H_Q3_3D_Outer[56] = {  0,  1,  2,  3,  4,  5,  6,  7,
                                    8,  9, 10, 11, 12, 13, 14, 15,
                                   16, 17, 18, 19, 20,         23,
                                   24,         27, 28, 29, 30, 31,
                                   32, 33, 34, 35, 36,         39,
                                   40,         43, 44, 45, 46, 47,
                                   48, 49, 50, 51, 52, 53, 54, 55,
                                   56, 57, 58, 59, 60, 61, 62, 63 };

static char C_H_Q3_3D_String[] = "C_H_Q3_3D";

TFEDesc3D *FE_C_H_Q3_3D_Obj=new TFEDesc3D(C_H_Q3_3D_String, C_H_Q3_3D_NDOF, 
                                C_H_Q3_3D_JointDOF,
                                C_H_Q3_3D_J, C_H_Q3_3D_NInner, C_H_Q3_3D_Inner,
                                C_H_Q3_3D_NOuter, C_H_Q3_3D_Outer);
