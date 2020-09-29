// hanging node for nonconforming P2/Q2 elements

static int N_P2_2D_N_Nodes_0 = 4;
static double N_P2_2D_Coupling_0[] = { 0.75, -0.25, -0.75, -0.25 };

THNDesc *HN_N_P2_2D_0_Obj = new THNDesc(N_P2_2D_N_Nodes_0, N_P2_2D_Coupling_0);

