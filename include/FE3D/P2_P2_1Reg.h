/*
    TFE3DMapper1Reg(char *name, char *description, int nfine, int ncoarse,
              int n_pairs, int **pairs,
              int n_hanging, int *hanging,
              HNDesc *hangingtypes, int **coupling,
              int n_nodes, int **twistpermutation);
*/
static char P2_P2_1Reg_Name[] = "P2_P2_1Reg";
static char P2_P2_1Reg_Desc[] = "conforming P2 element, 1-regular";
static int P2_P2_1Reg_NFine = 6;
static int P2_P2_1Reg_NCoarse = 6;
static int P2_P2_1Reg_N_Pairs = 15;
static int P2_P2_1Reg_Pairs0[][2] = { {0,24}, {2,11}, {4,21}, {5,14},
                                      {6,29}, {8,17}, {10,19}, {11,18},
                                      {12,26}, {14,23}, {16,22}, {17,20},
                                      {18,27}, {20,28}, {23,25} };
static int *P2_P2_1Reg_Pairs[1] = { (int *)P2_P2_1Reg_Pairs0 };

static int P2_P2_1Reg_NNodes = 30;

static int P2_P2_1Reg_NHanging = 9;
static int P2_P2_1Reg_Hanging[9] = { 1, 3, 4, 7, 9, 10, 13, 15, 16 };
static HNDesc P2_P2_1Reg_HangingTypes[9] = { HN_C_P2_3D_E, HN_C_P2_3D_E,
                                             HN_C_P2_3D_F, HN_C_P2_3D_E,
                                             HN_C_P2_3D_E, HN_C_P2_3D_F,
                                             HN_C_P2_3D_E, HN_C_P2_3D_E,
                                             HN_C_P2_3D_F };
static int P2_P2_1Reg_HN0[] = { 0, 2, 6 };
static int P2_P2_1Reg_HN1[] = { 0, 5, 12 };
static int P2_P2_1Reg_HN2[] = { 2, 6, 5, 8, 12 };
static int P2_P2_1Reg_HN3[] = { 6, 8, 12 };
static int P2_P2_1Reg_HN4[] = { 6, 2, 0 };
static int P2_P2_1Reg_HN5[] = { 8, 12, 2, 5, 0 };
static int P2_P2_1Reg_HN6[] = { 12, 5, 0 };
static int P2_P2_1Reg_HN7[] = { 12, 8, 6 };
static int P2_P2_1Reg_HN8[] = { 5, 0, 8, 2, 6 };

static int *P2_P2_1Reg_Coupling[9] = { P2_P2_1Reg_HN0, P2_P2_1Reg_HN1,
                                       P2_P2_1Reg_HN2, P2_P2_1Reg_HN3,
                                       P2_P2_1Reg_HN4, P2_P2_1Reg_HN5,
                                       P2_P2_1Reg_HN6, P2_P2_1Reg_HN7,
                                       P2_P2_1Reg_HN8 };

static int P2_P2_1Reg_TwistPerm0[] = { 0, 1, 2, 3, 4, 5 };
static int P2_P2_1Reg_TwistPerm1[] = { 2, 4, 5, 1, 3, 0 };
static int P2_P2_1Reg_TwistPerm2[] = { 5, 3, 0, 4, 1, 2 };

static int *P2_P2_1Reg_TwistPerm[3] = { P2_P2_1Reg_TwistPerm0,
                                        P2_P2_1Reg_TwistPerm1,
                                        P2_P2_1Reg_TwistPerm2 };

static int P2_P2_1Reg_NNoOpposite = 0;
static int **P2_P2_1Reg_NoOpposite = NULL;
        
TFE3DMapper1Reg *P2_P2_1Reg = new TFE3DMapper1Reg(
        P2_P2_1Reg_Name, P2_P2_1Reg_Desc,
        P2_P2_1Reg_NFine, P2_P2_1Reg_NCoarse,
        P2_P2_1Reg_N_Pairs, P2_P2_1Reg_Pairs,
        P2_P2_1Reg_NNoOpposite, P2_P2_1Reg_NoOpposite,
        P2_P2_1Reg_NHanging, P2_P2_1Reg_Hanging,
        P2_P2_1Reg_HangingTypes, P2_P2_1Reg_Coupling,
        P2_P2_1Reg_NNodes, P2_P2_1Reg_TwistPerm);
