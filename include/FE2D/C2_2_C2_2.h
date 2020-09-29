// P2 and Q2 like behaviour on considered edge
static char C2_2_C2_2_Name[] = "C2_2_C2_2";
static char C2_2_C2_2_Desc[] = "conforming P2 or Q2 element";
static int C2_2_C2_2_N0 = 3;
static int C2_2_C2_2_N1 = 3;
static int C2_2_C2_2_NPairs = 3;
static int C2_2_C2_2_Pairs[][2] = { {0,5}, {1,4}, {2,3} };
static int C2_2_C2_2_NHanging = 0;
static int *C2_2_C2_2_Hanging = NULL;
static HNDesc *C2_2_C2_2_HangingTypes = NULL;
static int **C2_2_C2_2_Coupling = NULL;
static int C2_2_C2_2_NFarHanging = 0;
static int *C2_2_C2_2_FarHanging = NULL;
static HNDesc *C2_2_C2_2_FarHangingTypes = NULL;
static int ****C2_2_C2_2_FarCoupling = NULL;
static int C2_2_C2_2_NNoopposite = 0;
static int *C2_2_C2_2_Nopposite = NULL;
static int C2_2_C2_2_NNodes = 6;

TFE2DMapper *C2_2_C2_2 = new TFE2DMapper(C2_2_C2_2_Name, C2_2_C2_2_Desc,
                             C2_2_C2_2_N0, C2_2_C2_2_N1,
                             C2_2_C2_2_NPairs, (int *)C2_2_C2_2_Pairs,
                             C2_2_C2_2_NHanging, C2_2_C2_2_Hanging,
                             C2_2_C2_2_HangingTypes, C2_2_C2_2_Coupling,
                             C2_2_C2_2_NFarHanging, C2_2_C2_2_FarHanging,
                             C2_2_C2_2_FarHangingTypes, C2_2_C2_2_FarCoupling,
                             C2_2_C2_2_NNoopposite, C2_2_C2_2_Nopposite,
                             C2_2_C2_2_NNodes);
