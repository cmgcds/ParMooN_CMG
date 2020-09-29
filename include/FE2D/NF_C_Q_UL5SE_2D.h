/*
TNodalFunctional2D(NodalFunctional2D id,
         int n_allfunctionals, int n_edgefunctionals, 
         int n_pointsall, int n_pointsedge, 
         double *xi, double *eta, double *t, 
         DoubleFunctVect *evalall, 
         DoubleFunctVect *evaledge); 
*/

static double NF_C_Q_UL5SE_2D_Xi[] = {
-1.0, -3.0/5.0, -1.0/5.0, 1.0/5.0, 3.0/5.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.0/5.0, 1.0/5.0,
-1.0/5.0, -3.0/5.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0/2.0, -1.0/2.0, -1.0/2.0, 0.0, 0.0, 0.0, 1.0/2.0, 1.0/2.0, 1.0/2.0,
};

static double NF_C_Q_UL5SE_2D_Eta[] = {
-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -3.0/5.0, -1.0/5.0, 1.0/5.0, 3.0/5.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 3.0/5.0, 1.0/5.0, -1.0/5.0, -3.0/5.0, -1.0/2.0, 0.0, 1.0/2.0, -1.0/2.0, 0.0, 1.0/2.0, -1.0/2.0, 0.0, 1.0/2.0,
};

static double NF_C_Q_UL5SE_2D_T[] = { -1.0, -3.0/5.0, -1.0/5.0, 1.0/5.0, 3.0/5.0, 1.0 };

void NF_C_Q_UL5SE_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  memcpy(Functionals, PointValues, 29*SizeOfDouble);
};

void NF_C_Q_UL5SE_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
  memcpy(Functionals, PointValues, 6*SizeOfDouble);
};

TNodalFunctional2D *NF_C_Q_UL5SE_2D_Obj = new TNodalFunctional2D
        (NF_C_Q_UL5SE_2D, 29, 6, 29, 6, NF_C_Q_UL5SE_2D_Xi, NF_C_Q_UL5SE_2D_Eta,
         NF_C_Q_UL5SE_2D_T, NF_C_Q_UL5SE_2D_EvalAll, NF_C_Q_UL5SE_2D_EvalEdge);
