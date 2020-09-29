/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

static double NF_C_T_P1_2D_Xi[] = { 0, 1, 0 };
static double NF_C_T_P1_2D_Eta[] = { 0, 0, 1 };
static double NF_C_T_P1_2D_T[] = { -1, 1 };

void NF_C_T_P1_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
}

void NF_C_T_P1_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
}

TNodalFunctional2D *NF_C_T_P1_2D_Obj = new TNodalFunctional2D
        (NF_C_T_P1_2D, 3, 2, 3, 2, NF_C_T_P1_2D_Xi, NF_C_T_P1_2D_Eta,
         NF_C_T_P1_2D_T, NF_C_T_P1_2D_EvalAll, NF_C_T_P1_2D_EvalEdge);
