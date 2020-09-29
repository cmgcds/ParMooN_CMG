

static double NF_C_L_P0_1D_Xi[] = { 0 };
static double NF_C_L_P0_1D_Eta[] = { 0 };
static double NF_C_L_P0_1D_T[] = { -1, 1 };

void NF_C_L_P0_1D_EvalAll(double *PointValues, double *Functionals)
{
  Functionals[0] = PointValues[0];
}

void NF_C_L_P0_1D_EvalEdge( double *PointValues, double *Functionals)
{

}

TNodalFunctional1D *NF_C_L_P0_1D_Obj = new TNodalFunctional1D
        (NF_C_L_P0_1D, 1, 0, 1, 0, NF_C_L_P0_1D_Xi, NF_C_L_P0_1D_Eta,
         NF_C_L_P0_1D_T, NF_C_L_P0_1D_EvalAll, NF_C_L_P0_1D_EvalEdge);

