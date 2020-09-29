// nonconfoming Q5 like behaviour on considered edge
static char N5_2_N5_2_Name[] = "N5_2_N5_2";
static char N5_2_N5_2_Desc[] = "nonconforming Q5 element";
static int N5_2_N5_2_N0 = 5;
static int N5_2_N5_2_N1 = 5;
static int N5_2_N5_2_NPairs = 5;
static int N5_2_N5_2_Pairs[][2] = { {0,5}, {1,6}, {2,7}, {3,8}, {4,9} };
static int N5_2_N5_2_NHanging = 0;
static int *N5_2_N5_2_Hanging = NULL;
static HNDesc *N5_2_N5_2_HangingTypes = NULL;
static int **N5_2_N5_2_Coupling = NULL;
static int N5_2_N5_2_NFarHanging = 0;
static int *N5_2_N5_2_FarHanging = NULL;
static HNDesc *N5_2_N5_2_FarHangingTypes = NULL;
static int ****N5_2_N5_2_FarCoupling = NULL;
static int N5_2_N5_2_NNoopposite = 0;
static int *N5_2_N5_2_Nopposite = NULL;
static int N5_2_N5_2_NNodes = 10;

TFE2DMapper *N5_2_N5_2 = new TFE2DMapper(N5_2_N5_2_Name, N5_2_N5_2_Desc,
                             N5_2_N5_2_N0, N5_2_N5_2_N1,
                             N5_2_N5_2_NPairs, (int *)N5_2_N5_2_Pairs,
                             N5_2_N5_2_NHanging, N5_2_N5_2_Hanging,
                             N5_2_N5_2_HangingTypes, N5_2_N5_2_Coupling,
                             N5_2_N5_2_NFarHanging, N5_2_N5_2_FarHanging,
                             N5_2_N5_2_FarHangingTypes, N5_2_N5_2_FarCoupling,
                             N5_2_N5_2_NNoopposite, N5_2_N5_2_Nopposite,
                             N5_2_N5_2_NNodes);
