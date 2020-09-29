// hanging node for conforming P1 elements in 3D

static int C_P1_3D_N_Nodes_E = 2;
static double C_P1_3D_Coupling_E[2] = { 0.5, 0.5 };

THNDesc *HN_C_P1_3D_E_Obj = new THNDesc(C_P1_3D_N_Nodes_E,
                                        C_P1_3D_Coupling_E);

