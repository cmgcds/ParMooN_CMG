// hanging node for conforming Q1 elements in 3D

static int C_Q1_3D_N_Nodes_E = 2;
static double C_Q1_3D_Coupling_E[2] = { 0.5, 0.5 };

THNDesc *HN_C_Q1_3D_E_Obj = new THNDesc(C_Q1_3D_N_Nodes_E,
                                        C_Q1_3D_Coupling_E);

static int C_Q1_3D_N_Nodes_F = 4;
static double C_Q1_3D_Coupling_F[4] = { 0.25, 0.25, 0.25, 0.25 };

THNDesc *HN_C_Q1_3D_F_Obj = new THNDesc(C_Q1_3D_N_Nodes_F,
                                        C_Q1_3D_Coupling_F);

