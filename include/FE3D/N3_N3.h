// mapper for hexahedron faces on both sides
static char N3_N3_Name[] = "N3_N3";
static char N3_N3_Desc[] = "nonconforming element of order 3";
static int N3_N3_N0 = 6;
static int N3_N3_N1 = 6;
static int N3_N3_NPairs = 6;
static int N3_N3_Pairs0[][2] = { {0,6}, {1,8}, {2,7}, {3,11}, {4,10}, {5,9} };
static int N3_N3_Pairs1[][2] = { {0,6}, {1,7}, {2,8}, {3,9}, {4,10}, {5,11} };
static int N3_N3_Pairs2[][2] = { {0,6}, {1,8}, {2,7}, {3,11}, {4,10}, {5,9} };
static int N3_N3_Pairs3[][2] = { {0,6}, {1,7}, {2,8}, {3,9}, {4,10}, {5,11} };
static int *N3_N3_Pairs[4] = { (int *)N3_N3_Pairs0, (int *)N3_N3_Pairs1,
                               (int *)N3_N3_Pairs2, (int *)N3_N3_Pairs3 };

static int N3_N3_NNodes = 12;

static int N3_N3_NNoOpposite = 0;
static int **N3_N3_NoOpposite = NULL;

TFE3DMapper *N3_N3 = new TFE3DMapper(N3_N3_Name, N3_N3_Desc,
                             N3_N3_N0, N3_N3_N1,
                             N3_N3_NPairs, N3_N3_Pairs,
                             N3_N3_NNoOpposite, N3_N3_NoOpposite,
                             N3_N3_NNodes);
