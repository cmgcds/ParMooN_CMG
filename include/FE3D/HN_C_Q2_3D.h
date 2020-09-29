// hanging node for conforming Q2 elements in 3D

static int C_Q2_3D_N_Nodes_E = 3;
static double C_Q2_3D_Coupling_E[3] = { 0.375, 0.75, -0.125 };

THNDesc *HN_C_Q2_3D_E_Obj = new THNDesc(C_Q2_3D_N_Nodes_E,
                                        C_Q2_3D_Coupling_E);

static int C_Q2_3D_N_Nodes_F = 9;
static double C_Q2_3D_Coupling_F[9] = {
        0.140625, 0.28125, -0.046875, 0.28125, 0.56250, -0.093750,
       -0.046875, -0.093750, 0.015625 };

THNDesc *HN_C_Q2_3D_F_Obj = new THNDesc(C_Q2_3D_N_Nodes_F,
                                        C_Q2_3D_Coupling_F);

