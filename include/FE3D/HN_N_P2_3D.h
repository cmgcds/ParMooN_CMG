// hanging node for conforming P2 elements in 3D

static int N_P2_3D_0_N_Nodes = 12;
static double N_P2_3D_0_Coupling[12] = { 0.25,  0.125, 0.125,
                                         0,     0,     0.125, 
                                         0,     0.125, 0,
                                         0.125, 0,     0.125 };

static int N_P2_3D_1_N_Nodes = 12;
static double N_P2_3D_1_Coupling[12] = { 0,     0,     0.125,
                                         0,     0.125, 0, 
                                         0.25,  0.125, 0.125,
                                         0,     0.125, 0.125 };

static int N_P2_3D_2_N_Nodes = 12;
static double N_P2_3D_2_Coupling[12] = { 0,     0.125, 0,
                                         0.25,  0.125, 0.125,
                                         0,     0,     0.125,
                                         0.125, 0.125, 0 };

THNDesc *HN_N_P2_3D_0_Obj = new THNDesc(N_P2_3D_0_N_Nodes,
                                        N_P2_3D_0_Coupling);

THNDesc *HN_N_P2_3D_1_Obj = new THNDesc(N_P2_3D_1_N_Nodes,
                                        N_P2_3D_1_Coupling);

THNDesc *HN_N_P2_3D_2_Obj = new THNDesc(N_P2_3D_2_N_Nodes,
                                        N_P2_3D_2_Coupling);

