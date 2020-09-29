// mapper for hexahedron faces on both sides
static char N1_N1_Name[] = "N1_N1";
static char N1_N1_Desc[] = "nonconforming Q1Rot element";
static int N1_N1_N0 = 1;
static int N1_N1_N1 = 1;
static int N1_N1_NPairs = 1;
static int N1_N1_Pairs0[][2] = { {0,1} };
static int N1_N1_Pairs1[][2] = { {0,1} };
static int N1_N1_Pairs2[][2] = { {0,1} };
static int N1_N1_Pairs3[][2] = { {0,1} };
static int *N1_N1_Pairs[4] = { (int *)N1_N1_Pairs0, (int *)N1_N1_Pairs1,
                               (int *)N1_N1_Pairs2, (int *)N1_N1_Pairs3 };

static int N1_N1_NNodes = 2;

static int N1_N1_NNoOpposite = 0;
static int **N1_N1_NoOpposite = NULL;

TFE3DMapper *N1_N1 = new TFE3DMapper(N1_N1_Name, N1_N1_Desc,
                             N1_N1_N0, N1_N1_N1,
                             N1_N1_NPairs, N1_N1_Pairs,
                             N1_N1_NNoOpposite, N1_N1_NoOpposite,
                             N1_N1_NNodes);
