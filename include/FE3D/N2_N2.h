// mapper for hexahedron faces on both sides
static char N2_N2_Name[] = "N2_N2";
static char N2_N2_Desc[] = "nonconforming element of order 2";
static int N2_N2_N0 = 3;
static int N2_N2_N1 = 3;
static int N2_N2_NPairs = 3;
static int N2_N2_Pairs0[][2] = { {0,3}, {1,5}, {2,4} };
static int N2_N2_Pairs1[][2] = { {0,3}, {1,4}, {2,5} };
static int N2_N2_Pairs2[][2] = { {0,3}, {1,5}, {2,4} };
static int N2_N2_Pairs3[][2] = { {0,3}, {1,4}, {2,5} };
static int *N2_N2_Pairs[4] = { (int *)N2_N2_Pairs0, (int *)N2_N2_Pairs1,
                               (int *)N2_N2_Pairs2, (int *)N2_N2_Pairs3 };

static int N2_N2_NNodes = 6;

static int N2_N2_NNoOpposite = 0;
static int **N2_N2_NoOpposite = NULL;

TFE3DMapper *N2_N2 = new TFE3DMapper(N2_N2_Name, N2_N2_Desc,
                             N2_N2_N0, N2_N2_N1,
                             N2_N2_NPairs, N2_N2_Pairs,
                             N2_N2_NNoOpposite, N2_N2_NoOpposite,
                             N2_N2_NNodes);
