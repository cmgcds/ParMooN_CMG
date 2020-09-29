// hanging node for conforming P4/Q4 elements

static int C_P4_2D_N_Nodes_0 = 5;
static double C_P4_2D_Coupling_0[] = {  -0.0390625, 0.21875, -0.546875, 1.09375, 0.2734375 };

static int C_P4_2D_N_Nodes_1 = 5;
static double C_P4_2D_Coupling_1[] = { 0.0234375, -0.15625, 0.703125, 0.46875, -0.0390625 };

static int C_P4_2D_N_Nodes_2 = 5;
static double C_P4_2D_Coupling_2[] = { -0.0390625, 0.46875, 0.703125, -0.15625, 0.0234375 };

static int C_P4_2D_N_Nodes_3 = 5;
static double C_P4_2D_Coupling_3[] = { 0.2734375, 1.09375, -0.546875, 0.21875, -0.0390625 };

THNDesc *HN_C_P4_2D_0_Obj = new THNDesc(C_P4_2D_N_Nodes_0, C_P4_2D_Coupling_0);
THNDesc *HN_C_P4_2D_1_Obj = new THNDesc(C_P4_2D_N_Nodes_1, C_P4_2D_Coupling_1);
THNDesc *HN_C_P4_2D_2_Obj = new THNDesc(C_P4_2D_N_Nodes_2, C_P4_2D_Coupling_2);
THNDesc *HN_C_P4_2D_3_Obj = new THNDesc(C_P4_2D_N_Nodes_3, C_P4_2D_Coupling_3);
