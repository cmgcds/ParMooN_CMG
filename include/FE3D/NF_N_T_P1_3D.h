/*
    TNodalFunctional3D(NodalFunctional3D id,
                       int n_allfunctionals, int *n_facefunctionals,
                       int B_T_Pointsall, int *B_T_Pointsface,
                       double *xi, double *eta, double *zeta,
                       double **xiarray, double **etaarray,
                       double **zetaarray,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evalface);
*/

/* for all functionals */
static double NF_N_T_P1_3D_Xi[]   = {  0.33333333333333333333,
                                     0.33333333333333333333,
                                     0.33333333333333333333,
                                     0 };

static double NF_N_T_P1_3D_Eta[]  = {  0.33333333333333333333,
                                     0,
                                     0.33333333333333333333,
                                     0.33333333333333333333 };
static double NF_N_T_P1_3D_Zeta[] = {  0,
                                     0.33333333333333333333,
                                     0.33333333333333333333,
                                     0.33333333333333333333 };

/* face 0                               0 */
static double NF_N_T_P1_3D_F0_Xi[]   = { 0.33333333333333333333 };
static double NF_N_T_P1_3D_F0_Eta[]  = { 0.33333333333333333333 };
static double NF_N_T_P1_3D_F0_Zeta[] = { 0 };

/* face 1                               1 */
static double NF_N_T_P1_3D_F1_Xi[]   = { 0.33333333333333333333 };
static double NF_N_T_P1_3D_F1_Eta[]  = { 0 };
static double NF_N_T_P1_3D_F1_Zeta[] = { 0.33333333333333333333 };

/* face 2                               2 */
static double NF_N_T_P1_3D_F2_Xi[]   = { 0.33333333333333333333 };
static double NF_N_T_P1_3D_F2_Eta[]  = { 0.33333333333333333333 };
static double NF_N_T_P1_3D_F2_Zeta[] = { 0.33333333333333333333 };

/* face 3                               3 */
static double NF_N_T_P1_3D_F3_Xi[]   = { 0 };
static double NF_N_T_P1_3D_F3_Eta[]  = { 0.33333333333333333333 };
static double NF_N_T_P1_3D_F3_Zeta[] = { 0.33333333333333333333 };

static double *NF_N_T_P1_3D_XiArray[6] = { 
                        NF_N_T_P1_3D_F0_Xi,
                        NF_N_T_P1_3D_F1_Xi,
                        NF_N_T_P1_3D_F2_Xi,
                        NF_N_T_P1_3D_F3_Xi };

static double *NF_N_T_P1_3D_EtaArray[4] = { 
                        NF_N_T_P1_3D_F0_Eta,
                        NF_N_T_P1_3D_F1_Eta,
                        NF_N_T_P1_3D_F2_Eta,
                        NF_N_T_P1_3D_F3_Eta };

static double *NF_N_T_P1_3D_ZetaArray[4] = { 
                        NF_N_T_P1_3D_F0_Zeta,
                        NF_N_T_P1_3D_F1_Zeta,
                        NF_N_T_P1_3D_F2_Zeta,
                        NF_N_T_P1_3D_F3_Zeta };

static double NF_N_T_P1_3D_T[1] = { 0.33333333333333333333 };
static double NF_N_T_P1_3D_S[1] = { 0.33333333333333333333 };

void NF_N_T_P1_3D_EvalAll(TCollection *Coll, TBaseCell *Cell,
                          double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 4*SizeOfDouble);
}

void NF_N_T_P1_3D_EvalFace(TCollection *Coll, TBaseCell *Cell, int Joint, 
                           double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
}

static int NF_N_T_P1_3D_N_AllFunctionals = 4;
static int NF_N_T_P1_3D_N_T_PointsAll = 4;
static int NF_N_T_P1_3D_N_FaceFunctionals[] = { 1, 1, 1, 1 };
static int NF_N_T_P1_3D_N_T_PointsFace[] = { 1, 1, 1, 1 };

TNodalFunctional3D *NF_N_T_P1_3D_Obj = new TNodalFunctional3D
        (NF_N_T_P1_3D, NF_N_T_P1_3D_N_AllFunctionals,
         NF_N_T_P1_3D_N_FaceFunctionals, NF_N_T_P1_3D_N_T_PointsAll,
         NF_N_T_P1_3D_N_T_PointsFace,
         NF_N_T_P1_3D_Xi, NF_N_T_P1_3D_Eta, NF_N_T_P1_3D_Zeta,
         NF_N_T_P1_3D_XiArray, NF_N_T_P1_3D_EtaArray,
         NF_N_T_P1_3D_ZetaArray,
         NF_N_T_P1_3D_T, NF_N_T_P1_3D_S,
         NF_N_T_P1_3D_EvalAll, NF_N_T_P1_3D_EvalFace);
