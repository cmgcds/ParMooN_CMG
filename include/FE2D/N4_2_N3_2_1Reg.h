// nonconforming P1 and Q1 like behaviour on considered edge
static char N4_2_N3_2_1Reg_Name[] = "N4_2_N3_2_1Reg";
static char N4_2_N3_2_1Reg_Desc[] = "nonconforming P4 or Q4 element, one regular grid";
static int N4_2_N3_2_1Reg_N0 = 4;
static int N4_2_N3_2_1Reg_N1 = 3;
static int N4_2_N3_2_1Reg_N2 = 3;
static int N4_2_N3_2_1Reg_NMid = 0;
static int *N4_2_N3_2_1Reg_Mid = NULL;
static int N4_2_N3_2_1Reg_NPairs = 0;
static int *N4_2_N3_2_1Reg_Pairs = NULL;
static int N4_2_N3_2_1Reg_NHanging = 3;
static int N4_2_N3_2_1Reg_Hanging[] = { 0, 1, 2 };
static HNDesc N4_2_N3_2_1Reg_HangingTypes[] = { HN_N_P1_2D_0, HN_N_P2_2D_0,
                                                HN_N_P3_2D_0 };
static int N4_2_N3_2_1Reg_Coupling_0[] = { 4,       7 };
static int N4_2_N3_2_1Reg_Coupling_1[] = { 4, 5,    7, 8 };
static int N4_2_N3_2_1Reg_Coupling_2[] = { 4, 5, 6, 7, 8, 9 };
static int *N4_2_N3_2_1Reg_Coupling[] = { N4_2_N3_2_1Reg_Coupling_0,
                                          N4_2_N3_2_1Reg_Coupling_1,
                                          N4_2_N3_2_1Reg_Coupling_2 };
static int N4_2_N3_2_1Reg_NFarHanging = 0;
static int *N4_2_N3_2_1Reg_FarHanging = NULL;
static HNDesc *N4_2_N3_2_1Reg_FarHangingTypes = NULL;
static int ****N4_2_N3_2_1Reg_FarCoupling = NULL;
static int N4_2_N3_2_1Reg_NNoopposite = 7;
static int N4_2_N3_2_1Reg_Nopposite[] = { 3, 4, 5, 6, 7, 8, 9 };
static int N4_2_N3_2_1Reg_NNodes = 10;

TFE2DMapper1Reg *N4_2_N3_2_1Reg = new TFE2DMapper1Reg(
                N4_2_N3_2_1Reg_Name, N4_2_N3_2_1Reg_Desc,
                N4_2_N3_2_1Reg_N0, N4_2_N3_2_1Reg_N1, N4_2_N3_2_1Reg_N2,
                N4_2_N3_2_1Reg_NPairs, (int *)N4_2_N3_2_1Reg_Pairs,
                N4_2_N3_2_1Reg_NMid, (int *)N4_2_N3_2_1Reg_Mid,
                N4_2_N3_2_1Reg_NHanging, N4_2_N3_2_1Reg_Hanging,
                N4_2_N3_2_1Reg_HangingTypes, N4_2_N3_2_1Reg_Coupling,
                N4_2_N3_2_1Reg_NFarHanging, N4_2_N3_2_1Reg_FarHanging,
                N4_2_N3_2_1Reg_FarHangingTypes, N4_2_N3_2_1Reg_FarCoupling,
                N4_2_N3_2_1Reg_NNoopposite, N4_2_N3_2_1Reg_Nopposite,
                N4_2_N3_2_1Reg_NNodes);
