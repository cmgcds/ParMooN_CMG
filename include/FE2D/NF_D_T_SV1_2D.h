/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

static double NF_D_T_SV1_2D_Xi[] = {
                                     7.0/18, 2.0/9, 13.0/18,
                                     7.0/18, 13.0/18, 2.0/9,
                                     2.0/9, 1.0/18, 1.0/18
                                   };
static double NF_D_T_SV1_2D_Eta[] = { 
                                      2.0/9, 1.0/18, 1.0/18,
                                      7.0/18, 2.0/9, 13.0/18,
                                      7.0/18, 13.0/18, 2.0/9
                                    };
static double *NF_D_T_SV1_2D_T = NULL;

void NF_D_T_SV1_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  Functionals[0] = ( 2*PointValues[0]+2*PointValues[1]+2*PointValues[2])/6;
  Functionals[1] = (-4*PointValues[0]-4*PointValues[1]+8*PointValues[2])/6;
  Functionals[2] = ( 8*PointValues[0]-4*PointValues[1]-4*PointValues[2])/6;
  
  Functionals[3] = ( 2*PointValues[3]+2*PointValues[4]+2*PointValues[5])/6;
  Functionals[4] = (-4*PointValues[3]-4*PointValues[4]+8*PointValues[5])/6;
  Functionals[5] = ( 8*PointValues[3]-4*PointValues[4]-4*PointValues[5])/6;

  Functionals[6] = ( 2*PointValues[6]+2*PointValues[7]+2*PointValues[8])/6;
  Functionals[7] = (-4*PointValues[6]-4*PointValues[7]+8*PointValues[8])/6;
  Functionals[8] = ( 8*PointValues[6]-4*PointValues[7]-4*PointValues[8])/6;
}

void NF_D_T_SV1_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
}

TNodalFunctional2D *NF_D_T_SV1_2D_Obj = new TNodalFunctional2D
        (NF_D_T_SV1_2D, 9, 0, 9, 0, NF_D_T_SV1_2D_Xi, NF_D_T_SV1_2D_Eta,
         NF_D_T_SV1_2D_T, NF_D_T_SV1_2D_EvalAll, NF_D_T_SV1_2D_EvalEdge);
