// ***********************************************************************
// P2 discontinuous element, 1D
// ***********************************************************************

// base function values
static void D_L_P2_1D_Funct(double xi, double *values)
{
  double t1 = xi*xi;
  double t2 = 2.5/3.;
  double t3 = 10./3.;
  double t4 = t2*xi;
  double t5 = 10.*t1;

  values[0]= 6.-15.*t1;
  values[1]= t3+t4-t5;
  values[2]= -t3+t4+t5;
}

// values of the derivatives in xi direction
static void D_L_P2_1D_DeriveXi(double xi, double *values)
{
 double t1 = 2.5/3.;
 double t2 = 20.*xi;

  values[0] = -30.*xi;
  values[1] = t1-t2;
  values[2] = t1+t2;
}

// values of the derivatives in xi-xi  direction
static void D_L_P2_1D_DeriveXiXi(double xi, double *values)
{
  values[0]=-30.;
  values[1]=-20.;
  values[2]=20.;
}

// ***********************************************************************
TBaseFunct1D *BF_D_L_P2_1D_Obj = new TBaseFunct1D
        (3, BF_D_L_P2_1D, D_L_P2_1D_Funct, D_L_P2_1D_DeriveXi, 
         D_L_P2_1D_DeriveXiXi, 2, 2);