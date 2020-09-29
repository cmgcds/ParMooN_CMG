// P5 and Q5 like behaviour on considered edge
static char C5_2_C5_2_Name[] = "C5_2_C5_2";
static char C5_2_C5_2_Desc[] = "conforming P5 or Q5 element";
static int C5_2_C5_2_N0 = 6;
static int C5_2_C5_2_N1 = 6;
static int C5_2_C5_2_NPairs = 6;
static int C5_2_C5_2_Pairs[][2] = { {0,11}, {1,10}, {2,9}, {3,8},
                                    {4,7}, {5,6}};
static int C5_2_C5_2_NHanging = 0;
static int *C5_2_C5_2_Hanging = NULL;
static HNDesc *C5_2_C5_2_HangingTypes = NULL;
static int **C5_2_C5_2_Coupling = NULL;
static int C5_2_C5_2_NFarHanging = 0;
static int *C5_2_C5_2_FarHanging = NULL;
static HNDesc *C5_2_C5_2_FarHangingTypes = NULL;
static int ****C5_2_C5_2_FarCoupling = NULL;
static int C5_2_C5_2_NNoopposite = 0;
static int *C5_2_C5_2_Nopposite = NULL;
static int C5_2_C5_2_NNodes = 12;

TFE2DMapper *C5_2_C5_2 = new TFE2DMapper(C5_2_C5_2_Name, C5_2_C5_2_Desc,
                             C5_2_C5_2_N0, C5_2_C5_2_N1,
                             C5_2_C5_2_NPairs, (int *)C5_2_C5_2_Pairs,
                             C5_2_C5_2_NHanging, C5_2_C5_2_Hanging,
                             C5_2_C5_2_HangingTypes, C5_2_C5_2_Coupling,
                             C5_2_C5_2_NFarHanging, C5_2_C5_2_FarHanging,
                             C5_2_C5_2_FarHangingTypes, C5_2_C5_2_FarCoupling,
                             C5_2_C5_2_NNoopposite, C5_2_C5_2_Nopposite,
                             C5_2_C5_2_NNodes);
