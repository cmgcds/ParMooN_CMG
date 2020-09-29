static double NF_C_Q_EL1_2D_Xi[] = { -1, 1, 1, -1,
 -.5773502691896257645091486, .5773502691896257645091486,
 -.5773502691896257645091486, .5773502691896257645091486  };

static double NF_C_Q_EL1_2D_Eta[] = { -1, -1, 1, 1,
  -.5773502691896257645091486, -.5773502691896257645091486,
   .5773502691896257645091486,  .5773502691896257645091486 };

static double NF_C_Q_EL1_2D_T[] = { -1, 1 };

void NF_C_Q_EL1_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = PointValues[3];
  Functionals[4] = PointValues[4]+PointValues[5]+PointValues[6]+PointValues[7];
}

void NF_C_Q_EL1_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
}

/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

TNodalFunctional2D *NF_C_Q_EL1_2D_Obj = new TNodalFunctional2D
        (NF_C_Q_EL1_2D, 5, 2, 8, 2, NF_C_Q_EL1_2D_Xi, NF_C_Q_EL1_2D_Eta,
         NF_C_Q_EL1_2D_T, NF_C_Q_EL1_2D_EvalAll, NF_C_Q_EL1_2D_EvalEdge);
