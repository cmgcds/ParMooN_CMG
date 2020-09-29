// P3 and Q3 like behaviour on considered edge
static char C3_2_C3_2_1Reg_Name[] = "C3_2_C3_2_1Reg";
static char C3_2_C3_2_1Reg_Desc[] = "conforming P3 or Q3 element, one regular grid";
static int C3_2_C3_2_1Reg_N0 = 4;
static int C3_2_C3_2_1Reg_N1 = 4;
static int C3_2_C3_2_1Reg_N2 = 4;
static int C3_2_C3_2_1Reg_NMid = 1;
static int C3_2_C3_2_1Reg_Mid[][2] = { {7,8} };
static int C3_2_C3_2_1Reg_NPairs = 4;
static int C3_2_C3_2_1Reg_Pairs[][2] = { {0,11}, {1,9}, {2,6}, {3,4} };
static int C3_2_C3_2_1Reg_NHanging = 3;
static int C3_2_C3_2_1Reg_Hanging[] = { 5, 7, 10 };
static HNDesc C3_2_C3_2_1Reg_HangingTypes[] = { HN_C_P3_2D_0, HN_C_P3_2D_1,
                                                 HN_C_P3_2D_2 };
static int C3_2_C3_2_1Reg_Coupling_0[] = { 0, 1, 2, 3 };
static int C3_2_C3_2_1Reg_Coupling_1[] = { 0, 1, 2, 3 };
static int C3_2_C3_2_1Reg_Coupling_2[] = { 0, 1, 2, 3 };
static int *C3_2_C3_2_1Reg_Coupling[] = { C3_2_C3_2_1Reg_Coupling_0,
                                          C3_2_C3_2_1Reg_Coupling_1,
                                          C3_2_C3_2_1Reg_Coupling_2 };
static int C3_2_C3_2_1Reg_NFarHanging = 0;
static int *C3_2_C3_2_1Reg_FarHanging = NULL;
static HNDesc *C3_2_C3_2_1Reg_FarHangingTypes = NULL;
static int ****C3_2_C3_2_1Reg_FarCoupling = NULL;
static int C3_2_C3_2_1Reg_NNoopposite = 0;
static int *C3_2_C3_2_1Reg_Nopposite = NULL;
static int C3_2_C3_2_1Reg_NNodes = 12;

TFE2DMapper1Reg *C3_2_C3_2_1Reg = new TFE2DMapper1Reg(
                C3_2_C3_2_1Reg_Name, C3_2_C3_2_1Reg_Desc,
                C3_2_C3_2_1Reg_N0, C3_2_C3_2_1Reg_N1, C3_2_C3_2_1Reg_N2,
                C3_2_C3_2_1Reg_NPairs, (int *)C3_2_C3_2_1Reg_Pairs,
                C3_2_C3_2_1Reg_NMid, (int *)C3_2_C3_2_1Reg_Mid,
                C3_2_C3_2_1Reg_NHanging, C3_2_C3_2_1Reg_Hanging,
                C3_2_C3_2_1Reg_HangingTypes, C3_2_C3_2_1Reg_Coupling,
                C3_2_C3_2_1Reg_NFarHanging, C3_2_C3_2_1Reg_FarHanging,
                C3_2_C3_2_1Reg_FarHangingTypes, C3_2_C3_2_1Reg_FarCoupling,
                C3_2_C3_2_1Reg_NNoopposite, C3_2_C3_2_1Reg_Nopposite,
                C3_2_C3_2_1Reg_NNodes);
