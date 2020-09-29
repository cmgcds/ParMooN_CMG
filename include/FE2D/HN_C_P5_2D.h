// hanging node for conforming P5/Q5 elements

static int C_P5_2D_N_Nodes_0 = 6;
static double C_P5_2D_Coupling_0[] = { (double)7/256, (double)-45/256, (double)63/128, (double)-105/128, (double)315/256, (double)63/256 };

static int C_P5_2D_N_Nodes_1 = 6;
static double C_P5_2D_Coupling_1[] = { (double)-3/256, (double)21/256, (double)-35/128, (double)105/128, (double)105/256, (double)-7/256 };

static int C_P5_2D_N_Nodes_2 = 6;
static double C_P5_2D_Coupling_2[] = { (double)3/256, (double)-25/256, (double)75/128, (double)75/128, (double)-25/256, (double)3/256 };

static int C_P5_2D_N_Nodes_3 = 6;
static double C_P5_2D_Coupling_3[] = { (double)-7/256, (double)105/256, (double)105/128, (double)-35/128, (double)21/256, (double)-3/256 };

static int C_P5_2D_N_Nodes_4 = 6;
static double C_P5_2D_Coupling_4[] = { (double)63/256, (double)315/256, (double)-105/128, (double)63/128, (double)-45/256, (double)7/256 };

THNDesc *HN_C_P5_2D_0_Obj = new THNDesc(C_P5_2D_N_Nodes_0, C_P5_2D_Coupling_0);
THNDesc *HN_C_P5_2D_1_Obj = new THNDesc(C_P5_2D_N_Nodes_1, C_P5_2D_Coupling_1);
THNDesc *HN_C_P5_2D_2_Obj = new THNDesc(C_P5_2D_N_Nodes_2, C_P5_2D_Coupling_2);
THNDesc *HN_C_P5_2D_3_Obj = new THNDesc(C_P5_2D_N_Nodes_3, C_P5_2D_Coupling_3);
THNDesc *HN_C_P5_2D_4_Obj = new THNDesc(C_P5_2D_N_Nodes_4, C_P5_2D_Coupling_4);
