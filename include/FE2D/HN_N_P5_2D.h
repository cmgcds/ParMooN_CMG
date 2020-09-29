// hanging node for nonconforming P5/Q5 elements

static int N_P5_2D_N_Nodes_0 = 10;
static double N_P5_2D_Coupling_0[] = { 0,  0.1875, 0.5625, -0.28125, 0.03125,
                                       0, -0.1875, 0.5625,  0.28125, 0.03125 };

THNDesc *HN_N_P5_2D_0_Obj = new THNDesc(N_P5_2D_N_Nodes_0, N_P5_2D_Coupling_0);

