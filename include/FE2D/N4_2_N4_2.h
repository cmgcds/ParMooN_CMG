// nonconfoming Q4 like behaviour on considered edge
static char N4_2_N4_2_Name[] = "N4_2_N4_2";
static char N4_2_N4_2_Desc[] = "nonconforming Q4 element";
static int N4_2_N4_2_N0 = 4;
static int N4_2_N4_2_N1 = 4;
static int N4_2_N4_2_NPairs = 4;
static int N4_2_N4_2_Pairs[][2] = { {0,4}, {1,5}, {2,6}, {3,7} };
static int N4_2_N4_2_NHanging = 0;
static int *N4_2_N4_2_Hanging = NULL;
static HNDesc *N4_2_N4_2_HangingTypes = NULL;
static int **N4_2_N4_2_Coupling = NULL;
static int N4_2_N4_2_NFarHanging = 0;
static int *N4_2_N4_2_FarHanging = NULL;
static HNDesc *N4_2_N4_2_FarHangingTypes = NULL;
static int ****N4_2_N4_2_FarCoupling = NULL;
static int N4_2_N4_2_NNoopposite = 0;
static int *N4_2_N4_2_Nopposite = NULL;
static int N4_2_N4_2_NNodes = 8;

TFE2DMapper *N4_2_N4_2 = new TFE2DMapper(N4_2_N4_2_Name, N4_2_N4_2_Desc,
                             N4_2_N4_2_N0, N4_2_N4_2_N1,
                             N4_2_N4_2_NPairs, (int *)N4_2_N4_2_Pairs,
                             N4_2_N4_2_NHanging, N4_2_N4_2_Hanging,
                             N4_2_N4_2_HangingTypes, N4_2_N4_2_Coupling,
                             N4_2_N4_2_NFarHanging, N4_2_N4_2_FarHanging,
                             N4_2_N4_2_FarHangingTypes, N4_2_N4_2_FarCoupling,
                             N4_2_N4_2_NNoopposite, N4_2_N4_2_Nopposite,
                             N4_2_N4_2_NNodes);
