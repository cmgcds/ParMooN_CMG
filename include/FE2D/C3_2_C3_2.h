// P3 and Q3 like behaviour on considered edge
static char C3_2_C3_2_Name[] = "C3_2_C3_2";
static char C3_2_C3_2_Desc[] = "conforming P3 or Q3 element";
static int C3_2_C3_2_N0 = 4;
static int C3_2_C3_2_N1 = 4;
static int C3_2_C3_2_NPairs = 4;
static int C3_2_C3_2_Pairs[][2] = { {0,7}, {1,6}, {2,5}, {3,4} };
static int C3_2_C3_2_NHanging = 0;
static int *C3_2_C3_2_Hanging = NULL;
static HNDesc *C3_2_C3_2_HangingTypes = NULL;
static int **C3_2_C3_2_Coupling = NULL;
static int C3_2_C3_2_NFarHanging = 0;
static int *C3_2_C3_2_FarHanging = NULL;
static HNDesc *C3_2_C3_2_FarHangingTypes = NULL;
static int ****C3_2_C3_2_FarCoupling = NULL;
static int C3_2_C3_2_NNoopposite = 0;
static int *C3_2_C3_2_Nopposite = NULL;
static int C3_2_C3_2_NNodes = 8;

TFE2DMapper *C3_2_C3_2 = new TFE2DMapper(C3_2_C3_2_Name, C3_2_C3_2_Desc,
                             C3_2_C3_2_N0, C3_2_C3_2_N1,
                             C3_2_C3_2_NPairs, (int *)C3_2_C3_2_Pairs,
                             C3_2_C3_2_NHanging, C3_2_C3_2_Hanging,
                             C3_2_C3_2_HangingTypes, C3_2_C3_2_Coupling,
                             C3_2_C3_2_NFarHanging, C3_2_C3_2_FarHanging,
                             C3_2_C3_2_FarHangingTypes, C3_2_C3_2_FarCoupling,
                             C3_2_C3_2_NNoopposite, C3_2_C3_2_Nopposite,
                             C3_2_C3_2_NNodes);
