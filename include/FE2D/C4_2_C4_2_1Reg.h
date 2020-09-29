// P4 and Q4 like behaviour on considered edge
static char C4_2_C4_2_1Reg_Name[] = "C4_2_C4_2_1Reg";
static char C4_2_C4_2_1Reg_Desc[] = "conforming P4 or Q4 element, one regular grid";
static int C4_2_C4_2_1Reg_N0 = 5;
static int C4_2_C4_2_1Reg_N1 = 5;
static int C4_2_C4_2_1Reg_N2 = 5;
static int C4_2_C4_2_1Reg_NMid = 1;
static int C4_2_C4_2_1Reg_Mid[][2] = { {9,10} };
static int C4_2_C4_2_1Reg_NPairs = 5;
static int C4_2_C4_2_1Reg_Pairs[][2] = { {0,14}, {1,12}, {2,9}, {3,7}, {4,5} };
static int C4_2_C4_2_1Reg_NHanging = 4;
static int C4_2_C4_2_1Reg_Hanging[] = { 6, 8, 11, 13 };
static HNDesc C4_2_C4_2_1Reg_HangingTypes[] = { HN_C_P4_2D_0, HN_C_P4_2D_1,
                                                 HN_C_P4_2D_2, HN_C_P4_2D_3 };
static int C4_2_C4_2_1Reg_Coupling_0[] = { 0, 1, 2, 3, 4 };
static int C4_2_C4_2_1Reg_Coupling_1[] = { 0, 1, 2, 3, 4 };
static int C4_2_C4_2_1Reg_Coupling_2[] = { 0, 1, 2, 3, 4 };
static int C4_2_C4_2_1Reg_Coupling_3[] = { 0, 1, 2, 3, 4 };
static int *C4_2_C4_2_1Reg_Coupling[] = { C4_2_C4_2_1Reg_Coupling_0,
                                          C4_2_C4_2_1Reg_Coupling_1,
                                          C4_2_C4_2_1Reg_Coupling_2,
      					  C4_2_C4_2_1Reg_Coupling_3 };
static int C4_2_C4_2_1Reg_NFarHanging = 0;
static int *C4_2_C4_2_1Reg_FarHanging = NULL;
static HNDesc *C4_2_C4_2_1Reg_FarHangingTypes = NULL;
static int ****C4_2_C4_2_1Reg_FarCoupling = NULL;
static int C4_2_C4_2_1Reg_NNoopposite = 0;
static int *C4_2_C4_2_1Reg_Nopposite = NULL;
static int C4_2_C4_2_1Reg_NNodes = 15;

TFE2DMapper1Reg *C4_2_C4_2_1Reg = new TFE2DMapper1Reg(
                C4_2_C4_2_1Reg_Name, C4_2_C4_2_1Reg_Desc,
                C4_2_C4_2_1Reg_N0, C4_2_C4_2_1Reg_N1, C4_2_C4_2_1Reg_N2,
                C4_2_C4_2_1Reg_NPairs, (int *)C4_2_C4_2_1Reg_Pairs,
                C4_2_C4_2_1Reg_NMid, (int *)C4_2_C4_2_1Reg_Mid,
                C4_2_C4_2_1Reg_NHanging, C4_2_C4_2_1Reg_Hanging,
                C4_2_C4_2_1Reg_HangingTypes, C4_2_C4_2_1Reg_Coupling,
                C4_2_C4_2_1Reg_NFarHanging, C4_2_C4_2_1Reg_FarHanging,
                C4_2_C4_2_1Reg_FarHangingTypes, C4_2_C4_2_1Reg_FarCoupling,
                C4_2_C4_2_1Reg_NNoopposite, C4_2_C4_2_1Reg_Nopposite,
                C4_2_C4_2_1Reg_NNodes);
