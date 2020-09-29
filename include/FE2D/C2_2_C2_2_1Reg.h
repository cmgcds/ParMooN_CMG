// P2 and Q2 like behaviour on considered edge
static char C2_2_C2_2_1Reg_Name[] = "C2_2_C2_2_1Reg";
static char C2_2_C2_2_1Reg_Desc[] = "conforming P2 or Q2 element, one regular grid";
static int C2_2_C2_2_1Reg_N0 = 3;
static int C2_2_C2_2_1Reg_N1 = 3;
static int C2_2_C2_2_1Reg_N2 = 3;
static int C2_2_C2_2_1Reg_NMid = 1;
static int C2_2_C2_2_1Reg_Mid[][2] = { {5,6} };
static int C2_2_C2_2_1Reg_NPairs = 3;
static int C2_2_C2_2_1Reg_Pairs[][2] = { {0,8}, {1,5}, {2,3} };
static int C2_2_C2_2_1Reg_NHanging = 2;
static int C2_2_C2_2_1Reg_Hanging[] = { 4, 7 };
static HNDesc C2_2_C2_2_1Reg_HangingTypes[] = { HN_C_P2_2D_0, HN_C_P2_2D_1 };
static int C2_2_C2_2_1Reg_Coupling_0[] = { 0, 1, 2 };
static int C2_2_C2_2_1Reg_Coupling_1[] = { 0, 1, 2 };
static int *C2_2_C2_2_1Reg_Coupling[] = { C2_2_C2_2_1Reg_Coupling_0,
                                          C2_2_C2_2_1Reg_Coupling_1 };
static int C2_2_C2_2_1Reg_NFarHanging = 0;
static int *C2_2_C2_2_1Reg_FarHanging = NULL;
static HNDesc *C2_2_C2_2_1Reg_FarHangingTypes = NULL;
static int ****C2_2_C2_2_1Reg_FarCoupling = NULL;
static int C2_2_C2_2_1Reg_NNoopposite = 0;
static int *C2_2_C2_2_1Reg_Nopposite = NULL;
static int C2_2_C2_2_1Reg_NNodes = 9;

TFE2DMapper1Reg *C2_2_C2_2_1Reg = new TFE2DMapper1Reg(
                C2_2_C2_2_1Reg_Name, C2_2_C2_2_1Reg_Desc,
                C2_2_C2_2_1Reg_N0, C2_2_C2_2_1Reg_N1, C2_2_C2_2_1Reg_N2,
                C2_2_C2_2_1Reg_NPairs, (int *)C2_2_C2_2_1Reg_Pairs,
                C2_2_C2_2_1Reg_NMid, (int *)C2_2_C2_2_1Reg_Mid,
                C2_2_C2_2_1Reg_NHanging, C2_2_C2_2_1Reg_Hanging,
                C2_2_C2_2_1Reg_HangingTypes, C2_2_C2_2_1Reg_Coupling,
                C2_2_C2_2_1Reg_NFarHanging, C2_2_C2_2_1Reg_FarHanging,
                C2_2_C2_2_1Reg_FarHangingTypes, C2_2_C2_2_1Reg_FarCoupling,
                C2_2_C2_2_1Reg_NNoopposite, C2_2_C2_2_1Reg_Nopposite,
                C2_2_C2_2_1Reg_NNodes);
