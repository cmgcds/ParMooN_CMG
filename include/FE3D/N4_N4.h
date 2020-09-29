// mapper for hexahedron faces on both sides
static char N4_N4_Name[] = "N4_N4";
static char N4_N4_Desc[] = "nonconforming element of order 4";
static int N4_N4_N0 = 10;
static int N4_N4_N1 = 10;
static int N4_N4_NPairs = 10;
static int N4_N4_Pairs0[][2] = { {0,10}, {1,12}, {2,11}, {3,15}, {4,14},
                                 {5,13}, {6,19}, {7,18}, {8,17}, {9,16} };
static int N4_N4_Pairs1[][2] = { {0,10}, {1,11}, {2,12}, {3,13}, {4,14},
                                 {5,15}, {6,16}, {7,17}, {8,18}, {9,19} };
static int N4_N4_Pairs2[][2] = { {0,10}, {1,12}, {2,11}, {3,15}, {4,14},
                                 {5,13}, {6,19}, {7,18}, {8,17}, {9,16} };
static int N4_N4_Pairs3[][2] = { {0,10}, {1,11}, {2,12}, {3,13}, {4,14},
                                 {5,15}, {6,16}, {7,17}, {8,18}, {9,19} };
static int *N4_N4_Pairs[4] = { (int *)N4_N4_Pairs0, (int *)N4_N4_Pairs1,
                               (int *)N4_N4_Pairs2, (int *)N4_N4_Pairs3 };

static int N4_N4_NNodes = 20;

static int N4_N4_NNoOpposite = 0;
static int **N4_N4_NoOpposite = NULL;

TFE3DMapper *N4_N4 = new TFE3DMapper(N4_N4_Name, N4_N4_Desc,
                             N4_N4_N0, N4_N4_N1,
                             N4_N4_NPairs, N4_N4_Pairs,
                             N4_N4_NNoOpposite, N4_N4_NoOpposite,
                             N4_N4_NNodes);
