static double NF_C_T_P1MINI_2D_Xi[7] = { 0, 0.5, 1, 0, 0.5, 0,
                                    0.33333333333333333333333 };
static double NF_C_T_P1MINI_2D_Eta[7] = { 0, 0, 0, 0.5, 0.5, 1,
                                     0.33333333333333333333333 };
static double NF_C_T_P1MINI_2D_T[2] = { -1, 1 };

void NF_C_T_P1MINI_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[2];
  Functionals[2] = PointValues[5];

  Functionals[3] = (27*PointValues[6]
                  + 8*(PointValues[1]+PointValues[3]+PointValues[4])
                  + 3*(PointValues[0]+PointValues[2]+PointValues[5]))/27;
}

void NF_C_T_P1MINI_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
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

TNodalFunctional2D *NF_C_T_P1MINI_2D_Obj = new TNodalFunctional2D
        (NF_C_T_P1MINI_2D, 4, 2, 7, 2,
         NF_C_T_P1MINI_2D_Xi, NF_C_T_P1MINI_2D_Eta,
         NF_C_T_P1MINI_2D_T, NF_C_T_P1MINI_2D_EvalAll,
         NF_C_T_P1MINI_2D_EvalEdge);

