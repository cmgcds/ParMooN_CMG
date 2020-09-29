// mapper for tetrahedron faces on both sides
static char NP1_NP1_Name[] = "NP1_NP1";
static char NP1_NP1_Desc[] = "nonconforming P1 element";
static int NP1_NP1_N0 = 1;
static int NP1_NP1_NP1 = 1;
static int NP1_NP1_NPairs = 1;
static int NP1_NP1_Pairs0[][2] = { {0,1} };
static int NP1_NP1_Pairs1[][2] = { {0,1} };
static int NP1_NP1_Pairs2[][2] = { {0,1} };
static int *NP1_NP1_Pairs[3] = { (int *)NP1_NP1_Pairs0,
                                 (int *)NP1_NP1_Pairs1,
                                 (int *)NP1_NP1_Pairs2 };

static int NP1_NP1_NNodes = 2;

static int NP1_NP1_NNoOpposite = 0;
static int **NP1_NP1_NoOpposite = NULL;

TFE3DMapper *NP1_NP1 = new TFE3DMapper(NP1_NP1_Name, NP1_NP1_Desc,
                             NP1_NP1_N0, NP1_NP1_NP1,
                             NP1_NP1_NPairs, NP1_NP1_Pairs,
                             NP1_NP1_NNoOpposite, NP1_NP1_NoOpposite,
                             NP1_NP1_NNodes);
