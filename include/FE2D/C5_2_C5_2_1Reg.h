// P5 and Q5 like behaviour on considered edge
static char C5_2_C5_2_1Reg_Name[] = "C5_2_C5_2_1Reg";
static char C5_2_C5_2_1Reg_Desc[] = "conforming P5 or Q5 element, one regular grid";
static int C5_2_C5_2_1Reg_N0 = 6;
static int C5_2_C5_2_1Reg_N1 = 6;
static int C5_2_C5_2_1Reg_N2 = 6;
static int C5_2_C5_2_1Reg_NMid = 1;
static int C5_2_C5_2_1Reg_Mid[][2] = { {11,12} };
static int C5_2_C5_2_1Reg_NPairs = 6;
static int C5_2_C5_2_1Reg_Pairs[][2] = { {0,17}, {1,15}, {2,13}, {3,10}, {4,8}, {5,6} };
static int C5_2_C5_2_1Reg_NHanging = 5;
static int C5_2_C5_2_1Reg_Hanging[] = { 7, 9, 11, 14, 16 };
static HNDesc C5_2_C5_2_1Reg_HangingTypes[] = { HN_C_P5_2D_0, HN_C_P5_2D_1,
                                                 HN_C_P5_2D_2, HN_C_P5_2D_3, HN_C_P5_2D_4 };
static int C5_2_C5_2_1Reg_Coupling_0[] = { 0, 1, 2, 3, 4, 5 };
static int C5_2_C5_2_1Reg_Coupling_1[] = { 0, 1, 2, 3, 4, 5 };
static int C5_2_C5_2_1Reg_Coupling_2[] = { 0, 1, 2, 3, 4, 5 };
static int C5_2_C5_2_1Reg_Coupling_3[] = { 0, 1, 2, 3, 4, 5 };
static int C5_2_C5_2_1Reg_Coupling_4[] = { 0, 1, 2, 3, 4, 5 };
static int *C5_2_C5_2_1Reg_Coupling[] = { C5_2_C5_2_1Reg_Coupling_0,
                                          C5_2_C5_2_1Reg_Coupling_1,
                                          C5_2_C5_2_1Reg_Coupling_2,
					  C5_2_C5_2_1Reg_Coupling_3,
					  C5_2_C5_2_1Reg_Coupling_4 };
static int C5_2_C5_2_1Reg_NFarHanging = 0;
static int *C5_2_C5_2_1Reg_FarHanging = NULL;
static HNDesc *C5_2_C5_2_1Reg_FarHangingTypes = NULL;
static int ****C5_2_C5_2_1Reg_FarCoupling = NULL;
static int C5_2_C5_2_1Reg_NNoopposite = 0;
static int *C5_2_C5_2_1Reg_Nopposite = NULL;
static int C5_2_C5_2_1Reg_NNodes = 18;

TFE2DMapper1Reg *C5_2_C5_2_1Reg = new TFE2DMapper1Reg(
                C5_2_C5_2_1Reg_Name, C5_2_C5_2_1Reg_Desc,
                C5_2_C5_2_1Reg_N0, C5_2_C5_2_1Reg_N1, C5_2_C5_2_1Reg_N2,
                C5_2_C5_2_1Reg_NPairs, (int *)C5_2_C5_2_1Reg_Pairs,
                C5_2_C5_2_1Reg_NMid, (int *)C5_2_C5_2_1Reg_Mid,
                C5_2_C5_2_1Reg_NHanging, C5_2_C5_2_1Reg_Hanging,
                C5_2_C5_2_1Reg_HangingTypes, C5_2_C5_2_1Reg_Coupling,
                C5_2_C5_2_1Reg_NFarHanging, C5_2_C5_2_1Reg_FarHanging,
                C5_2_C5_2_1Reg_FarHangingTypes, C5_2_C5_2_1Reg_FarCoupling,
                C5_2_C5_2_1Reg_NNoopposite, C5_2_C5_2_1Reg_Nopposite,
                C5_2_C5_2_1Reg_NNodes);
