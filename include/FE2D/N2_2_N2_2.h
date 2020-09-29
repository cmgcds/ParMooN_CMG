// nonconfoming Q2 like behaviour on considered edge
static char N2_2_N2_2_Name[] = "N2_2_N2_2";
static char N2_2_N2_2_Desc[] = "nonconforming Q2 element";
static int N2_2_N2_2_N0 = 2;
static int N2_2_N2_2_N1 = 2;
static int N2_2_N2_2_NPairs = 2;
static int N2_2_N2_2_Pairs[][2] = { {0,2}, {1,3} };
static int N2_2_N2_2_NHanging = 0;
static int *N2_2_N2_2_Hanging = NULL;
static HNDesc *N2_2_N2_2_HangingTypes = NULL;
static int **N2_2_N2_2_Coupling = NULL;
static int N2_2_N2_2_NFarHanging = 0;
static int *N2_2_N2_2_FarHanging = NULL;
static HNDesc *N2_2_N2_2_FarHangingTypes = NULL;
static int ****N2_2_N2_2_FarCoupling = NULL;
static int N2_2_N2_2_NNoopposite = 0;
static int *N2_2_N2_2_Nopposite = NULL;
static int N2_2_N2_2_NNodes = 4;

TFE2DMapper *N2_2_N2_2 = new TFE2DMapper(N2_2_N2_2_Name, N2_2_N2_2_Desc,
                             N2_2_N2_2_N0, N2_2_N2_2_N1,
                             N2_2_N2_2_NPairs, (int *)N2_2_N2_2_Pairs,
                             N2_2_N2_2_NHanging, N2_2_N2_2_Hanging,
                             N2_2_N2_2_HangingTypes, N2_2_N2_2_Coupling,
                             N2_2_N2_2_NFarHanging, N2_2_N2_2_FarHanging,
                             N2_2_N2_2_FarHangingTypes, N2_2_N2_2_FarCoupling,
                             N2_2_N2_2_NNoopposite, N2_2_N2_2_Nopposite,
                             N2_2_N2_2_NNodes);
