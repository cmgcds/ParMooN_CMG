/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

static double NF_C_Q_Q0_2D_Xi[] = { 0 };
static double NF_C_Q_Q0_2D_Eta[] = { 0 };
static double *NF_C_Q_Q0_2D_T = NULL;

void NF_C_Q_Q0_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  Functionals[0] = PointValues[0];
}

void NF_C_Q_Q0_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
}

TNodalFunctional2D *NF_C_Q_Q0_2D_Obj = new TNodalFunctional2D
        (NF_C_Q_Q0_2D, 1, 0, 1, 0, NF_C_Q_Q0_2D_Xi, NF_C_Q_Q0_2D_Eta,
         NF_C_Q_Q0_2D_T, NF_C_Q_Q0_2D_EvalAll, NF_C_Q_Q0_2D_EvalEdge);
