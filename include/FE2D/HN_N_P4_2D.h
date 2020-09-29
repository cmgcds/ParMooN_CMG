// hanging node for nonconforming P4/Q4 elements

static int N_P4_2D_N_Nodes_0 = 8;
static double N_P4_2D_Coupling_0[] = { -0.4375, -0.4375,  0.4375, -0.0625,
                                        0.4375, -0.4375, -0.4375, -0.0625 };

THNDesc *HN_N_P4_2D_0_Obj = new THNDesc(N_P4_2D_N_Nodes_0, N_P4_2D_Coupling_0);

