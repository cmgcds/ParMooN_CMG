// mapper for tetrahedron faces on both sides
static char P2_P2_Name[] = "P2_P2";
static char P2_P2_Desc[] = "conforming P2 element";
static int P2_P2_N0 = 6;
static int P2_P2_N1 = 6;
static int P2_P2_NPairs = 6;
static int P2_P2_Pairs0[][2] = { {0,6}, {1,9}, {2,11},
                                 {3,7}, {4,10}, {5,8} };
static int P2_P2_Pairs1[][2] = { {0,8}, {1,7}, {2,6},
                                 {3,10}, {4,9}, {5,11} };
static int P2_P2_Pairs2[][2] = { {0,11}, {1,10}, {2,8},
                                 {3,9}, {4,7}, {5,6} };
static int *P2_P2_Pairs[3] = { (int *)P2_P2_Pairs0, (int *)P2_P2_Pairs1,
                               (int *)P2_P2_Pairs2 };

static int P2_P2_NNodes = 12;

static int P2_P2_NNoOpposite = 0;
static int **P2_P2_NoOpposite = NULL;

TFE3DMapper *P2_P2 = new TFE3DMapper(P2_P2_Name, P2_P2_Desc,
                             P2_P2_N0, P2_P2_N1,
                             P2_P2_NPairs, P2_P2_Pairs,
                             P2_P2_NNoOpposite, P2_P2_NoOpposite,
                             P2_P2_NNodes);
