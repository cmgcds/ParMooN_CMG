// hanging node for conforming P2/Q2 elements

static int C_P2_2D_N_Nodes_0 = 3;
static double C_P2_2D_Coupling_0[] = { -0.125, 0.75, 0.375 };

THNDesc *HN_C_P2_2D_0_Obj = new THNDesc(C_P2_2D_N_Nodes_0, C_P2_2D_Coupling_0);

static int C_P2_2D_N_Nodes_1 = 3;
static double C_P2_2D_Coupling_1[] = { 0.375, 0.75, -0.125 };

THNDesc *HN_C_P2_2D_1_Obj = new THNDesc(C_P2_2D_N_Nodes_1, C_P2_2D_Coupling_1);

