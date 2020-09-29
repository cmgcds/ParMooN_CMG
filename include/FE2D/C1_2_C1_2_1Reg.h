// P1 and Q1 like behaviour on considered edge
static char C1_2_C1_2_1Reg_Name[] = "C1_2_C1_2_1Reg";
static char C1_2_C1_2_1Reg_Desc[] = "conforming P1 or Q1 element, one regular grid";
static int C1_2_C1_2_1Reg_N0 = 2;
static int C1_2_C1_2_1Reg_N1 = 2;
static int C1_2_C1_2_1Reg_N2 = 2;
static int C1_2_C1_2_1Reg_NMid = 1;
static int C1_2_C1_2_1Reg_Mid[][2] = { {3,4} };
static int C1_2_C1_2_1Reg_NPairs = 2;
static int C1_2_C1_2_1Reg_Pairs[][2] = { {0,5}, {1,2} };
static int C1_2_C1_2_1Reg_NHanging = 1;
static int C1_2_C1_2_1Reg_Hanging[] = { 3 };
static HNDesc C1_2_C1_2_1Reg_HangingTypes[] = { HN_C_P1_2D_0 };
static int C1_2_C1_2_1Reg_Coupling_0[] = { 0, 1 };
static int *C1_2_C1_2_1Reg_Coupling[] = { C1_2_C1_2_1Reg_Coupling_0 };
static int C1_2_C1_2_1Reg_NFarHanging = 0;
static int *C1_2_C1_2_1Reg_FarHanging = NULL;
static HNDesc *C1_2_C1_2_1Reg_FarHangingTypes = NULL;
static int ****C1_2_C1_2_1Reg_FarCoupling = NULL;
static int C1_2_C1_2_1Reg_NNoopposite = 0;
static int *C1_2_C1_2_1Reg_Nopposite = NULL;
static int C1_2_C1_2_1Reg_NNodes = 6;

TFE2DMapper1Reg *C1_2_C1_2_1Reg = new TFE2DMapper1Reg(
                C1_2_C1_2_1Reg_Name, C1_2_C1_2_1Reg_Desc,
                C1_2_C1_2_1Reg_N0, C1_2_C1_2_1Reg_N1, C1_2_C1_2_1Reg_N2,
                C1_2_C1_2_1Reg_NPairs, (int *)C1_2_C1_2_1Reg_Pairs,
                C1_2_C1_2_1Reg_NMid, (int *)C1_2_C1_2_1Reg_Mid,
                C1_2_C1_2_1Reg_NHanging, C1_2_C1_2_1Reg_Hanging,
                C1_2_C1_2_1Reg_HangingTypes, C1_2_C1_2_1Reg_Coupling,
                C1_2_C1_2_1Reg_NFarHanging, C1_2_C1_2_1Reg_FarHanging,
                C1_2_C1_2_1Reg_FarHangingTypes, C1_2_C1_2_1Reg_FarCoupling,
                C1_2_C1_2_1Reg_NNoopposite, C1_2_C1_2_1Reg_Nopposite,
                C1_2_C1_2_1Reg_NNodes);
