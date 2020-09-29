/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

static double NF_D_T_P1_2D_Xi[] = { 0.5, 0.5, 0 };
static double NF_D_T_P1_2D_Eta[] = { 0, 0.5, 0.5 };
static double NF_D_T_P1_2D_T_P[] = 
        { -0.77459666924148337703585307995647992, 0,
           0.77459666924148337703585307995647992 };

/*
   weighting functions: 2, 24*xi-8, 24*eta-8
*/

void NF_D_T_P1_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  Functionals[0] = 2*(PointValues[0]+PointValues[1]+PointValues[2])/6;
  Functionals[1] = ( 4*PointValues[0]
                    +4*PointValues[1]
                    -8*PointValues[2])/6;
  Functionals[2] = (-8*PointValues[0]
                    +4*PointValues[1]
                    +4*PointValues[2])/6;
}

void NF_D_T_P1_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
}

TNodalFunctional2D *NF_D_T_P1_2D_Obj = new TNodalFunctional2D
        (NF_D_T_P1_2D, 3, 0, 3, 0, NF_D_T_P1_2D_Xi, NF_D_T_P1_2D_Eta,
         NF_D_T_P1_2D_T_P, NF_D_T_P1_2D_EvalAll, NF_D_T_P1_2D_EvalEdge);

