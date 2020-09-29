// P0 and Q0 like behaviour on considered edge
static char C0_2_C0_2_Name[] = "C0_2_C0_2";
static char C0_2_C0_2_Desc[] = "conforming P0 or Q0 element";
static int C0_2_C0_2_N0 = 0;
static int C0_2_C0_2_N1 = 0;
static int C0_2_C0_2_NPairs = 0;
static int *C0_2_C0_2_Pairs = NULL;
static int C0_2_C0_2_NHanging = 0;
static int *C0_2_C0_2_Hanging = NULL;
static HNDesc *C0_2_C0_2_HangingTypes = NULL;
static int **C0_2_C0_2_Coupling = NULL;
static int C0_2_C0_2_NFarHanging = 0;
static int *C0_2_C0_2_FarHanging = NULL;
static HNDesc *C0_2_C0_2_FarHangingTypes = NULL;
static int ****C0_2_C0_2_FarCoupling = NULL;
static int C0_2_C0_2_NNoopposite = 0;
static int *C0_2_C0_2_Nopposite = NULL;
static int C0_2_C0_2_NNodes = 0;

TFE2DMapper *C0_2_C0_2 = new TFE2DMapper(C0_2_C0_2_Name, C0_2_C0_2_Desc,
                             C0_2_C0_2_N0, C0_2_C0_2_N1,
                             C0_2_C0_2_NPairs, (int *)C0_2_C0_2_Pairs,
                             C0_2_C0_2_NHanging, C0_2_C0_2_Hanging,
                             C0_2_C0_2_HangingTypes, C0_2_C0_2_Coupling,
                             C0_2_C0_2_NFarHanging, C0_2_C0_2_FarHanging,
                             C0_2_C0_2_FarHangingTypes, C0_2_C0_2_FarCoupling,
                             C0_2_C0_2_NNoopposite, C0_2_C0_2_Nopposite,
                             C0_2_C0_2_NNodes);
