// hanging node for conforming P2 elements in 3D

static int C_P2_3D_N_Nodes_E = 3;
static double C_P2_3D_Coupling_E[3] = { 0.375, 0.75, -0.125 };

THNDesc *HN_C_P2_3D_E_Obj = new THNDesc(C_P2_3D_N_Nodes_E,
                                        C_P2_3D_Coupling_E);

static int C_P2_3D_N_Nodes_F = 5;
static double C_P2_3D_Coupling_F[5] = { 0.5, -0.125, 0.5, 0.25, -0.125 };

THNDesc *HN_C_P2_3D_F_Obj = new THNDesc(C_P2_3D_N_Nodes_F,
                                        C_P2_3D_Coupling_F);

