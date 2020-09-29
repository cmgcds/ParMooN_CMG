// P0 and Q0 like behaviour on considered edge
static char C0_2_C0_2_1Reg_Name[] = "C0_2_C0_2_1Reg";
static char C0_2_C0_2_1Reg_Desc[] = "conforming P0 or Q0 element, one regular grid";
static int C0_2_C0_2_1Reg_N0 = 0;
static int C0_2_C0_2_1Reg_N1 = 0;
static int C0_2_C0_2_1Reg_N2 = 0;
static int C0_2_C0_2_1Reg_NMid = 0;
static int *C0_2_C0_2_1Reg_Mid = NULL;
static int C0_2_C0_2_1Reg_NPairs = 0;
static int *C0_2_C0_2_1Reg_Pairs = NULL;
static int C0_2_C0_2_1Reg_NHanging = 0;
static int *C0_2_C0_2_1Reg_Hanging = NULL;
static HNDesc *C0_2_C0_2_1Reg_HangingTypes = NULL;
static int **C0_2_C0_2_1Reg_Coupling = NULL;
static int C0_2_C0_2_1Reg_NFarHanging = 0;
static int *C0_2_C0_2_1Reg_FarHanging = NULL;
static HNDesc *C0_2_C0_2_1Reg_FarHangingTypes = NULL;
static int ****C0_2_C0_2_1Reg_FarCoupling = NULL;
static int C0_2_C0_2_1Reg_NNoopposite = 0;
static int *C0_2_C0_2_1Reg_Nopposite = NULL;
static int C0_2_C0_2_1Reg_NNodes = 0;

TFE2DMapper1Reg *C0_2_C0_2_1Reg = new TFE2DMapper1Reg(
                C0_2_C0_2_1Reg_Name, C0_2_C0_2_1Reg_Desc,
                C0_2_C0_2_1Reg_N0, C0_2_C0_2_1Reg_N1, C0_2_C0_2_1Reg_N2,
                C0_2_C0_2_1Reg_NPairs, (int *)C0_2_C0_2_1Reg_Pairs,
                C0_2_C0_2_1Reg_NMid, (int *)C0_2_C0_2_1Reg_Mid,
                C0_2_C0_2_1Reg_NHanging, C0_2_C0_2_1Reg_Hanging,
                C0_2_C0_2_1Reg_HangingTypes, C0_2_C0_2_1Reg_Coupling,
                C0_2_C0_2_1Reg_NFarHanging, C0_2_C0_2_1Reg_FarHanging,
                C0_2_C0_2_1Reg_FarHangingTypes, C0_2_C0_2_1Reg_FarCoupling,
                C0_2_C0_2_1Reg_NNoopposite, C0_2_C0_2_1Reg_Nopposite,
                C0_2_C0_2_1Reg_NNodes);
