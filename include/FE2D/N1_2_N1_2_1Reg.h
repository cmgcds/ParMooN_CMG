// nonconforming P1 and Q1 like behaviour on considered edge
static char N1_2_N1_2_1Reg_Name[] = "N1_2_N1_2_1Reg";
static char N1_2_N1_2_1Reg_Desc[] = "nonconforming P1 or Q1 element, one regular grid";
static int N1_2_N1_2_1Reg_N0 = 1;
static int N1_2_N1_2_1Reg_N1 = 1;
static int N1_2_N1_2_1Reg_N2 = 1;
static int N1_2_N1_2_1Reg_NMid = 0;
static int *N1_2_N1_2_1Reg_Mid = NULL;
static int N1_2_N1_2_1Reg_NPairs = 0;
static int *N1_2_N1_2_1Reg_Pairs = NULL;
static int N1_2_N1_2_1Reg_NHanging = 1;
static int N1_2_N1_2_1Reg_Hanging[] = { 0 };
static HNDesc N1_2_N1_2_1Reg_HangingTypes[] = { HN_N_P1_2D_0 };
static int N1_2_N1_2_1Reg_Coupling_0[] = { 1, 2 };
static int *N1_2_N1_2_1Reg_Coupling[] = { N1_2_N1_2_1Reg_Coupling_0 };
static int N1_2_N1_2_1Reg_NFarHanging = 0;
static int *N1_2_N1_2_1Reg_FarHanging = NULL;
static HNDesc *N1_2_N1_2_1Reg_FarHangingTypes = NULL;
static int ****N1_2_N1_2_1Reg_FarCoupling = NULL;
static int N1_2_N1_2_1Reg_NNoopposite = 2;
static int N1_2_N1_2_1Reg_Nopposite[] = { 1, 2 };
static int N1_2_N1_2_1Reg_NNodes = 3;

TFE2DMapper1Reg *N1_2_N1_2_1Reg = new TFE2DMapper1Reg(
                N1_2_N1_2_1Reg_Name, N1_2_N1_2_1Reg_Desc,
                N1_2_N1_2_1Reg_N0, N1_2_N1_2_1Reg_N1, N1_2_N1_2_1Reg_N2,
                N1_2_N1_2_1Reg_NPairs, (int *)N1_2_N1_2_1Reg_Pairs,
                N1_2_N1_2_1Reg_NMid, (int *)N1_2_N1_2_1Reg_Mid,
                N1_2_N1_2_1Reg_NHanging, N1_2_N1_2_1Reg_Hanging,
                N1_2_N1_2_1Reg_HangingTypes, N1_2_N1_2_1Reg_Coupling,
                N1_2_N1_2_1Reg_NFarHanging, N1_2_N1_2_1Reg_FarHanging,
                N1_2_N1_2_1Reg_FarHangingTypes, N1_2_N1_2_1Reg_FarCoupling,
                N1_2_N1_2_1Reg_NNoopposite, N1_2_N1_2_1Reg_Nopposite,
                N1_2_N1_2_1Reg_NNodes);
