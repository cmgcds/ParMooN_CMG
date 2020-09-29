/*
    TFE3DMapper1Reg(char *name, char *description, int nfine, int ncoarse,
              int n_pairs, int **pairs,
              int n_hanging, int *hanging,
              HNDesc *hangingtypes, int **coupling,
              int n_nodes, int **twistpermutation);
*/
static char NP1_NP1_1Reg_Name[] = "NP1_NP1_1Reg";
static char NP1_NP1_1Reg_Desc[] = "conforming NP1 element, 1-regular";
static int NP1_NP1_1Reg_NFine = 1;
static int NP1_NP1_1Reg_NCoarse = 1;
static int NP1_NP1_1Reg_N_Pairs = 0;
static int *NP1_NP1_1Reg_Pairs0 = NULL;
static int *NP1_NP1_1Reg_Pairs[1] = { (int *)NP1_NP1_1Reg_Pairs0 };

static int NP1_NP1_1Reg_NNodes = 5;

static int NP1_NP1_1Reg_NHanging = 1;
static int NP1_NP1_1Reg_Hanging[1] = { 4 };
static HNDesc NP1_NP1_1Reg_HangingTypes[1] = { HN_N_P1_3D_E };

static int NP1_NP1_1Reg_HN0[] = { 0, 1, 2, 3 };

static int *NP1_NP1_1Reg_Coupling[1] = { NP1_NP1_1Reg_HN0 };

static int NP1_NP1_1Reg_TwistPerm0[] = { 0 };
static int NP1_NP1_1Reg_TwistPerm1[] = { 0 };
static int NP1_NP1_1Reg_TwistPerm2[] = { 0 };

static int *NP1_NP1_1Reg_TwistPerm[3] = { NP1_NP1_1Reg_TwistPerm0,
                                   NP1_NP1_1Reg_TwistPerm1,
                                   NP1_NP1_1Reg_TwistPerm2 };

static int NP1_NP1_1Reg_NNoOpposite = 4;
static int NP1_NP1_1Reg_NNoOpposite1[] = { 0, 1, 2, 3 };
static int *NP1_NP1_1Reg_NoOpposite[] = { NP1_NP1_1Reg_NNoOpposite1 };
        
TFE3DMapper1Reg *NP1_NP1_1Reg = new TFE3DMapper1Reg(
        NP1_NP1_1Reg_Name, NP1_NP1_1Reg_Desc,
        NP1_NP1_1Reg_NFine, NP1_NP1_1Reg_NCoarse,
        NP1_NP1_1Reg_N_Pairs, NP1_NP1_1Reg_Pairs,
        NP1_NP1_1Reg_NNoOpposite, NP1_NP1_1Reg_NoOpposite,
        NP1_NP1_1Reg_NHanging, NP1_NP1_1Reg_Hanging,
        NP1_NP1_1Reg_HangingTypes, NP1_NP1_1Reg_Coupling,
        NP1_NP1_1Reg_NNodes, NP1_NP1_1Reg_TwistPerm);
