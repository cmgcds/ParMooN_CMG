// nonconfoming Q3 like behaviour on considered edge
static char N3_2_N3_2_Name[] = "N3_2_N3_2";
static char N3_2_N3_2_Desc[] = "nonconforming Q3 element";
static int N3_2_N3_2_N0 = 3;
static int N3_2_N3_2_N1 = 3;
static int N3_2_N3_2_NPairs = 3;
static int N3_2_N3_2_Pairs[][2] = { {0,3}, {1,4}, {2,5} };
static int N3_2_N3_2_NHanging = 0;
static int *N3_2_N3_2_Hanging = NULL;
static HNDesc *N3_2_N3_2_HangingTypes = NULL;
static int **N3_2_N3_2_Coupling = NULL;
static int N3_2_N3_2_NFarHanging = 0;
static int *N3_2_N3_2_FarHanging = NULL;
static HNDesc *N3_2_N3_2_FarHangingTypes = NULL;
static int ****N3_2_N3_2_FarCoupling = NULL;
static int N3_2_N3_2_NNoopposite = 0;
static int *N3_2_N3_2_Nopposite = NULL;
static int N3_2_N3_2_NNodes = 6;

TFE2DMapper *N3_2_N3_2 = new TFE2DMapper(N3_2_N3_2_Name, N3_2_N3_2_Desc,
                             N3_2_N3_2_N0, N3_2_N3_2_N1,
                             N3_2_N3_2_NPairs, (int *)N3_2_N3_2_Pairs,
                             N3_2_N3_2_NHanging, N3_2_N3_2_Hanging,
                             N3_2_N3_2_HangingTypes, N3_2_N3_2_Coupling,
                             N3_2_N3_2_NFarHanging, N3_2_N3_2_FarHanging,
                             N3_2_N3_2_FarHangingTypes, N3_2_N3_2_FarCoupling,
                             N3_2_N3_2_NNoopposite, N3_2_N3_2_Nopposite,
                             N3_2_N3_2_NNodes);
