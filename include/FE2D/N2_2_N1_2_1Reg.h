// nonconforming P1 and Q1 like behaviour on considered edge
static char N2_2_N1_2_1Reg_Name[] = "N2_2_N1_2_1Reg";
static char N2_2_N1_2_1Reg_Desc[] = "nonconforming P2 or Q2 element, one regular grid";
static int N2_2_N1_2_1Reg_N0 = 2;
static int N2_2_N1_2_1Reg_N1 = 1;
static int N2_2_N1_2_1Reg_N2 = 1;
static int N2_2_N1_2_1Reg_NMid = 0;
static int *N2_2_N1_2_1Reg_Mid = NULL;
static int N2_2_N1_2_1Reg_NPairs = 0;
static int *N2_2_N1_2_1Reg_Pairs = NULL;
static int N2_2_N1_2_1Reg_NHanging = 1;
static int N2_2_N1_2_1Reg_Hanging[] = { 0 };
static HNDesc N2_2_N1_2_1Reg_HangingTypes[] = { HN_N_P1_2D_0 };
static int N2_2_N1_2_1Reg_Coupling_0[] = { 2, 3 };
static int *N2_2_N1_2_1Reg_Coupling[] = { N2_2_N1_2_1Reg_Coupling_0 };
static int N2_2_N1_2_1Reg_NFarHanging = 0;
static int *N2_2_N1_2_1Reg_FarHanging = NULL;
static HNDesc *N2_2_N1_2_1Reg_FarHangingTypes = NULL;
static int ****N2_2_N1_2_1Reg_FarCoupling = NULL;
static int N2_2_N1_2_1Reg_NNoopposite = 3;
static int N2_2_N1_2_1Reg_Nopposite[] = { 1, 2, 3};
static int N2_2_N1_2_1Reg_NNodes = 4;

TFE2DMapper1Reg *N2_2_N1_2_1Reg = new TFE2DMapper1Reg(
                N2_2_N1_2_1Reg_Name, N2_2_N1_2_1Reg_Desc,
                N2_2_N1_2_1Reg_N0, N2_2_N1_2_1Reg_N2, N2_2_N1_2_1Reg_N2,
                N2_2_N1_2_1Reg_NPairs, (int *)N2_2_N1_2_1Reg_Pairs,
                N2_2_N1_2_1Reg_NMid, (int *)N2_2_N1_2_1Reg_Mid,
                N2_2_N1_2_1Reg_NHanging, N2_2_N1_2_1Reg_Hanging,
                N2_2_N1_2_1Reg_HangingTypes, N2_2_N1_2_1Reg_Coupling,
                N2_2_N1_2_1Reg_NFarHanging, N2_2_N1_2_1Reg_FarHanging,
                N2_2_N1_2_1Reg_FarHangingTypes, N2_2_N1_2_1Reg_FarCoupling,
                N2_2_N1_2_1Reg_NNoopposite, N2_2_N1_2_1Reg_Nopposite,
                N2_2_N1_2_1Reg_NNodes);
