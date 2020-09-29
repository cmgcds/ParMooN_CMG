// mapper for hexahedron faces on both sides
static char Q2_Q2_Name[] = "Q2_Q2";
static char Q2_Q2_Desc[] = "conforming Q2 element";
static int Q2_Q2_N0 = 9;
static int Q2_Q2_N1 = 9;
static int Q2_Q2_NPairs = 9;
static int Q2_Q2_Pairs0[][2] = { {0,9}, {1,12}, {2,15}, {3,10}, {4,13}, {5,16},
                                 {6,11}, {7,14}, {8,17} };
static int Q2_Q2_Pairs1[][2] = { {0,11}, {1,10}, {2,9}, {3,14}, {4,13}, {5,12},
                                 {6,17}, {7,16}, {8,15} }; 
static int Q2_Q2_Pairs2[][2] = { {0,17}, {1,14}, {2,11}, {3,16}, {4,13}, {5,10},
                                 {6,15}, {7,12}, {8, 9} };
static int Q2_Q2_Pairs3[][2] = { {0,15}, {1,16}, {2,17}, {3,12}, {4,13}, {5,14},
                                 {6,9}, {7,10}, {8,11} };
static int *Q2_Q2_Pairs[4] = { (int *)Q2_Q2_Pairs0, (int *)Q2_Q2_Pairs1,
                               (int *)Q2_Q2_Pairs2, (int *)Q2_Q2_Pairs3 };

static int Q2_Q2_NNodes = 18;

static int Q2_Q2_NNoOpposite = 0;
static int **Q2_Q2_NoOpposite = NULL;

TFE3DMapper *Q2_Q2 = new TFE3DMapper(Q2_Q2_Name, Q2_Q2_Desc,
                             Q2_Q2_N0, Q2_Q2_N1,
                             Q2_Q2_NPairs, Q2_Q2_Pairs,
                             Q2_Q2_NNoOpposite, Q2_Q2_NoOpposite,
                             Q2_Q2_NNodes);
