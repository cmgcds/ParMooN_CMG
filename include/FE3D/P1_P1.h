// mapper for tetrahedron faces on both sides
static char P1_P1_Name[] = "P1_P1";
static char P1_P1_Desc[] = "conforming P1 element";
static int P1_P1_N0 = 3;
static int P1_P1_N1 = 3;
static int P1_P1_NPairs = 3;
static int P1_P1_Pairs0[][2] = { {0,3}, {1,5}, {2,4} };
static int P1_P1_Pairs1[][2] = { {0,4}, {1,3}, {2,5} };
static int P1_P1_Pairs2[][2] = { {0,5}, {1,4}, {2,3} };
static int *P1_P1_Pairs[3] = { (int *)P1_P1_Pairs0, (int *)P1_P1_Pairs1,
                               (int *)P1_P1_Pairs2 };

static int P1_P1_NNodes = 6;

static int P1_P1_NNoOpposite = 0;
static int **P1_P1_NoOpposite = NULL;

TFE3DMapper *P1_P1 = new TFE3DMapper(P1_P1_Name, P1_P1_Desc,
                             P1_P1_N0, P1_P1_N1,
                             P1_P1_NPairs, P1_P1_Pairs,
                             P1_P1_NNoOpposite, P1_P1_NoOpposite,
                             P1_P1_NNodes);
