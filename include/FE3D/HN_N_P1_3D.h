// hanging node for conforming P1 elements in 3D

static int N_P1_3D_N_Nodes_E = 4;
static double N_P1_3D_Coupling_E[4] = { 0.25, 0.25, 0.25, 0.25 };

THNDesc *HN_N_P1_3D_E_Obj = new THNDesc(N_P1_3D_N_Nodes_E,
                                        N_P1_3D_Coupling_E);

