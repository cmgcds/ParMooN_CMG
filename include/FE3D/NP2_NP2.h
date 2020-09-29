// mapper for tetrahedral faces on both sides
static char NP2_NP2_Name[] = "NP2_NP2";
static char NP2_NP2_Desc[] = "nonconforming P2 element";
static int NP2_NP2_N0 = 3;
static int NP2_NP2_NP2 = 3;
static int NP2_NP2_NPairs = 3;
static int NP2_NP2_Pairs0[][2] = { {0,3}, {1,5}, {2,4} };
static int NP2_NP2_Pairs1[][2] = { {0,4}, {1,3}, {2,5} };
static int NP2_NP2_Pairs2[][2] = { {0,5}, {1,4}, {2,3} };
static int *NP2_NP2_Pairs[3] = { (int *)NP2_NP2_Pairs0,
                                 (int *)NP2_NP2_Pairs1,
                                 (int *)NP2_NP2_Pairs2 };

static int NP2_NP2_NNodes = 6;

static int NP2_NP2_NNoOpposite = 0;
static int **NP2_NP2_NoOpposite = NULL;

TFE3DMapper *NP2_NP2 = new TFE3DMapper(NP2_NP2_Name, NP2_NP2_Desc,
                             NP2_NP2_N0, NP2_NP2_NP2,
                             NP2_NP2_NPairs, NP2_NP2_Pairs,
                             NP2_NP2_NNoOpposite, NP2_NP2_NoOpposite,
                             NP2_NP2_NNodes);
