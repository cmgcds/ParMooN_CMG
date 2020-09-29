// P1 and Q1 like behaviour on considered edge
static char C1_2_C1_2_Name[] = "C1_2_C1_2";
static char C1_2_C1_2_Desc[] = "conforming P1 or Q1 element";
static int C1_2_C1_2_N0 = 2;
static int C1_2_C1_2_N1 = 2;
static int C1_2_C1_2_NPairs = 2;
static int C1_2_C1_2_Pairs[][2] = { {0,3}, {1,2} };
static int C1_2_C1_2_NHanging = 0;
static int *C1_2_C1_2_Hanging = NULL;
static HNDesc *C1_2_C1_2_HangingTypes = NULL;
static int **C1_2_C1_2_Coupling = NULL;
static int C1_2_C1_2_NFarHanging = 0;
static int *C1_2_C1_2_FarHanging = NULL;
static HNDesc *C1_2_C1_2_FarHangingTypes = NULL;
static int ****C1_2_C1_2_FarCoupling = NULL;
static int C1_2_C1_2_NNoopposite = 0;
static int *C1_2_C1_2_Nopposite = NULL;
static int C1_2_C1_2_NNodes = 4;

TFE2DMapper *C1_2_C1_2 = new TFE2DMapper(C1_2_C1_2_Name, C1_2_C1_2_Desc,
                             C1_2_C1_2_N0, C1_2_C1_2_N1,
                             C1_2_C1_2_NPairs, (int *)C1_2_C1_2_Pairs,
                             C1_2_C1_2_NHanging, C1_2_C1_2_Hanging,
                             C1_2_C1_2_HangingTypes, C1_2_C1_2_Coupling,
                             C1_2_C1_2_NFarHanging, C1_2_C1_2_FarHanging,
                             C1_2_C1_2_FarHangingTypes, C1_2_C1_2_FarCoupling,
                             C1_2_C1_2_NNoopposite, C1_2_C1_2_Nopposite,
                             C1_2_C1_2_NNodes);
