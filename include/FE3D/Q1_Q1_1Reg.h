/*
    TFE3DMapper1Reg(char *name, char *description, int nfine, int ncoarse,
              int n_pairs, int **pairs,
              int n_hanging, int *hanging,
              HNDesc *hangingtypes, int **coupling,
              int n_nodes, int **twistpermutation);
*/
static char Q1_Q1_1Reg_Name[] = "Q1_Q1_1Reg";
static char Q1_Q1_1Reg_Desc[] = "conforming Q1 element, 1-regular";
static int Q1_Q1_1Reg_NFine = 4;
static int Q1_Q1_1Reg_NCoarse = 4;
static int Q1_Q1_1Reg_N_Pairs = 11;
static int Q1_Q1_1Reg_Pairs0[][2] = { {0,16}, { 1,6}, {2,13}, {3,7},
                                      {4,18}, {5,10}, {7,11}, {8,19},
                                      {9,14}, {11,15}, {12,17} };
static int *Q1_Q1_1Reg_Pairs[1] = { (int *)Q1_Q1_1Reg_Pairs0 };

static int Q1_Q1_1Reg_NNodes = 20;

static int Q1_Q1_1Reg_NHanging = 5;
static int Q1_Q1_1Reg_Hanging[5] = { 1, 5, 9, 2, 3 };
static HNDesc Q1_Q1_1Reg_HangingTypes[5] = { HN_C_Q1_3D_E, HN_C_Q1_3D_E,
                                             HN_C_Q1_3D_E, HN_C_Q1_3D_E,
                                             HN_C_Q1_3D_F };
static int Q1_Q1_1Reg_HN0[] = { 0, 4 };
static int Q1_Q1_1Reg_HN1[] = { 4, 8 };
static int Q1_Q1_1Reg_HN2[] = { 8, 12 };
static int Q1_Q1_1Reg_HN3[] = { 0, 12 };
static int Q1_Q1_1Reg_HN4[] = { 0, 4, 8, 12 };

static int *Q1_Q1_1Reg_Coupling[5] = { Q1_Q1_1Reg_HN0, Q1_Q1_1Reg_HN1,
                                       Q1_Q1_1Reg_HN2, Q1_Q1_1Reg_HN3,
                                       Q1_Q1_1Reg_HN4 };

static int Q1_Q1_1Reg_TwistPerm0[] = { 0, 1, 2, 3 };
static int Q1_Q1_1Reg_TwistPerm1[] = { 1, 3, 0, 2 };
static int Q1_Q1_1Reg_TwistPerm2[] = { 3, 2, 1, 0 };
static int Q1_Q1_1Reg_TwistPerm3[] = { 2, 0, 3, 1 };

static int *Q1_Q1_1Reg_TwistPerm[4] = { Q1_Q1_1Reg_TwistPerm0,
                                   Q1_Q1_1Reg_TwistPerm1,
                                   Q1_Q1_1Reg_TwistPerm2,
                                   Q1_Q1_1Reg_TwistPerm3 };

static int Q1_Q1_1Reg_NNoOpposite = 0;
static int **Q1_Q1_1Reg_NoOpposite = NULL;

TFE3DMapper1Reg *Q1_Q1_1Reg = new TFE3DMapper1Reg(
        Q1_Q1_1Reg_Name, Q1_Q1_1Reg_Desc,
        Q1_Q1_1Reg_NFine, Q1_Q1_1Reg_NCoarse,
        Q1_Q1_1Reg_N_Pairs, Q1_Q1_1Reg_Pairs,
        Q1_Q1_1Reg_NNoOpposite, Q1_Q1_1Reg_NoOpposite,
        Q1_Q1_1Reg_NHanging, Q1_Q1_1Reg_Hanging,
        Q1_Q1_1Reg_HangingTypes, Q1_Q1_1Reg_Coupling,
        Q1_Q1_1Reg_NNodes, Q1_Q1_1Reg_TwistPerm);
