/*
    TNodalFunctional3D(NodalFunctional3D id,
                       int n_allfunctionals, int *n_facefunctionals,
                       int n_pointsall, int *n_pointsface,
                       double *xi, double *eta, double *zeta,
                       double **xiarray, double **etaarray,
                       double **zetaarray,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evalface);
*/

/* for all functionals */
static double NF_C_H_Q00_3D_Xi[]   = {  0 };
static double NF_C_H_Q00_3D_Eta[]  = {  0 };
static double NF_C_H_Q00_3D_Zeta[] = {  0 };

/* face 0                               0 */
static double *NF_C_H_Q00_3D_F0_Xi = NULL;
static double *NF_C_H_Q00_3D_F0_Eta = NULL;
static double *NF_C_H_Q00_3D_F0_Zeta = NULL;

/* face 1                               1 */
static double *NF_C_H_Q00_3D_F1_Xi = NULL;
static double *NF_C_H_Q00_3D_F1_Eta = NULL;
static double *NF_C_H_Q00_3D_F1_Zeta = NULL;

/* face 2                               2 */
static double *NF_C_H_Q00_3D_F2_Xi = NULL;
static double *NF_C_H_Q00_3D_F2_Eta = NULL;
static double *NF_C_H_Q00_3D_F2_Zeta = NULL;

/* face 3                               3 */
static double *NF_C_H_Q00_3D_F3_Xi = NULL;
static double *NF_C_H_Q00_3D_F3_Eta = NULL;
static double *NF_C_H_Q00_3D_F3_Zeta = NULL;

/* face 4                               4 */
static double *NF_C_H_Q00_3D_F4_Xi = NULL;
static double *NF_C_H_Q00_3D_F4_Eta = NULL;
static double *NF_C_H_Q00_3D_F4_Zeta = NULL;

/* face 5                               5 */
static double *NF_C_H_Q00_3D_F5_Xi = NULL;
static double *NF_C_H_Q00_3D_F5_Eta = NULL;
static double *NF_C_H_Q00_3D_F5_Zeta = NULL;

static double *NF_C_H_Q00_3D_XiArray[6] = { 
                        NF_C_H_Q00_3D_F0_Xi,
                        NF_C_H_Q00_3D_F1_Xi,
                        NF_C_H_Q00_3D_F2_Xi,
                        NF_C_H_Q00_3D_F3_Xi,
                        NF_C_H_Q00_3D_F4_Xi,
                        NF_C_H_Q00_3D_F5_Xi };

static double *NF_C_H_Q00_3D_EtaArray[6] = { 
                        NF_C_H_Q00_3D_F0_Eta,
                        NF_C_H_Q00_3D_F1_Eta,
                        NF_C_H_Q00_3D_F2_Eta,
                        NF_C_H_Q00_3D_F3_Eta,
                        NF_C_H_Q00_3D_F4_Eta,
                        NF_C_H_Q00_3D_F5_Eta };

static double *NF_C_H_Q00_3D_ZetaArray[6] = { 
                        NF_C_H_Q00_3D_F0_Zeta,
                        NF_C_H_Q00_3D_F1_Zeta,
                        NF_C_H_Q00_3D_F2_Zeta,
                        NF_C_H_Q00_3D_F3_Zeta,
                        NF_C_H_Q00_3D_F4_Zeta,
                        NF_C_H_Q00_3D_F5_Zeta };

static double *NF_C_H_Q00_3D_T = NULL;
static double *NF_C_H_Q00_3D_S = NULL;

void NF_C_H_Q00_3D_EvalAll(TCollection *Coll, TBaseCell *Cell,
                           double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
}

void NF_C_H_Q00_3D_EvalFace(TCollection *Coll, TBaseCell *Cell, int Joint, 
                            double *PointValues, double *Functionals)
{
}

static int NF_C_H_Q00_3D_N_AllFunctionals = 1;
static int NF_C_H_Q00_3D_N_PointsAll = 1;
static int NF_C_H_Q00_3D_N_FaceFunctionals[] = { 0, 0, 0, 0, 0, 0 };
static int NF_C_H_Q00_3D_N_PointsFace[] = { 0, 0, 0, 0, 0, 0 };

TNodalFunctional3D *NF_C_H_Q00_3D_Obj = new TNodalFunctional3D
        (NF_C_H_Q00_3D, NF_C_H_Q00_3D_N_AllFunctionals,
         NF_C_H_Q00_3D_N_FaceFunctionals, NF_C_H_Q00_3D_N_PointsAll,
         NF_C_H_Q00_3D_N_PointsFace,
         NF_C_H_Q00_3D_Xi, NF_C_H_Q00_3D_Eta, NF_C_H_Q00_3D_Zeta,
         NF_C_H_Q00_3D_XiArray, NF_C_H_Q00_3D_EtaArray,
         NF_C_H_Q00_3D_ZetaArray,
         NF_C_H_Q00_3D_T, NF_C_H_Q00_3D_S,
         NF_C_H_Q00_3D_EvalAll, NF_C_H_Q00_3D_EvalFace);
