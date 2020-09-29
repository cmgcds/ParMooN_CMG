// P9 and Q9 like behaviour on considered edge
static char C9_2_C9_2_Name[] = "C9_2_C9_2";
static char C9_2_C9_2_Desc[] = "conforming P9 or Q9 element";
static int C9_2_C9_2_N0 = 10;
static int C9_2_C9_2_N1 = 10;
static int C9_2_C9_2_NPairs = 10;
static int C9_2_C9_2_Pairs[][2] = { {0,19}, {1,18}, {2,17}, {3,16}, {4,15},
                                    {5,14}, {6,13}, {7,12}, {8,11}, {9,10} };
static int C9_2_C9_2_NHanging = 0;
static int *C9_2_C9_2_Hanging = NULL;
static HNDesc *C9_2_C9_2_HangingTypes = NULL;
static int **C9_2_C9_2_Coupling = NULL;
static int C9_2_C9_2_NFarHanging = 0;
static int *C9_2_C9_2_FarHanging = NULL;
static HNDesc *C9_2_C9_2_FarHangingTypes = NULL;
static int ****C9_2_C9_2_FarCoupling = NULL;
static int C9_2_C9_2_NNoopposite = 0;
static int *C9_2_C9_2_Nopposite = NULL;
static int C9_2_C9_2_NNodes = 20;

TFE2DMapper *C9_2_C9_2 = new TFE2DMapper(C9_2_C9_2_Name, C9_2_C9_2_Desc,
                             C9_2_C9_2_N0, C9_2_C9_2_N1,
                             C9_2_C9_2_NPairs, (int *)C9_2_C9_2_Pairs,
                             C9_2_C9_2_NHanging, C9_2_C9_2_Hanging,
                             C9_2_C9_2_HangingTypes, C9_2_C9_2_Coupling,
                             C9_2_C9_2_NFarHanging, C9_2_C9_2_FarHanging,
                             C9_2_C9_2_FarHangingTypes, C9_2_C9_2_FarCoupling,
                             C9_2_C9_2_NNoopposite, C9_2_C9_2_Nopposite,
                             C9_2_C9_2_NNodes);
