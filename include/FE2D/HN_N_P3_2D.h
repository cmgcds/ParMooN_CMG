// hanging node for nonconforming P3/Q3 elements

static int N_P3_2D_N_Nodes_0 = 6;
static double N_P3_2D_Coupling_0[] = { 0, -0.625, 0.125,
                                       0,  0.625, 0.125 };

THNDesc *HN_N_P3_2D_0_Obj = new THNDesc(N_P3_2D_N_Nodes_0, N_P3_2D_Coupling_0);

