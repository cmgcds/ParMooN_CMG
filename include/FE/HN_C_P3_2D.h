// hanging node for conforming P3/Q3 elements

static int C_P3_2D_N_Nodes_0 = 4;
static double C_P3_2D_Coupling_0[] = {  0.0625, -0.3125,  0.9375,  0.3125 };

static int C_P3_2D_N_Nodes_1 = 4;
static double C_P3_2D_Coupling_1[] = { -0.0625,  0.5625,  0.5625, -0.0625 };

static int C_P3_2D_N_Nodes_2 = 4;
static double C_P3_2D_Coupling_2[] = {  0.3125,  0.9375, -0.3125,  0.0625 };

THNDesc *HN_C_P3_2D_0_Obj = new THNDesc(C_P3_2D_N_Nodes_0, C_P3_2D_Coupling_0);
THNDesc *HN_C_P3_2D_1_Obj = new THNDesc(C_P3_2D_N_Nodes_1, C_P3_2D_Coupling_1);
THNDesc *HN_C_P3_2D_2_Obj = new THNDesc(C_P3_2D_N_Nodes_2, C_P3_2D_Coupling_2);

