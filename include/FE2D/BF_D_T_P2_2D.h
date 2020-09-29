// ***********************************************************************
// P2 element, discontinous, 2D, triangle
// ***********************************************************************

// base function values
static void D_T_P2_2D_Funct(double xi, double eta, double *values)
{
  double t7, t9, t10, t11;

  t7 = xi*xi;
  t9 = xi*eta;
  t10 = 6.0*t9;
  t11 = eta*eta;

  values[0] = 1.0;
  values[1] = -8.0+24.0*xi;
  values[2] = -8.0+24.0*eta;
  values[3] = 1.0-6.0*xi-2.0*eta+6.0*t7+t10+t11;
  values[4] = 1.0-4.0*xi-4.0*eta+3.0*t7+8.0*t9+3.0*t11;
  values[5] = 1.0-2.0*xi-6.0*eta+t7+t10+6.0*t11;
}

// values of the derivatives in xi direction
static void D_T_P2_2D_DeriveXi(double xi, double eta, double *values)
{
  double t2;

  t2 = 6.0*eta;

  values[0] = 0.0;
  values[1] = 24.0;
  values[2] = 0.0;
  values[3] = -6.0+12.0*xi+t2;
  values[4] = -4.0+6.0*xi+8.0*eta;
  values[5] = -2.0+2.0*xi+t2;
}

// values of the derivatives in eta direction
static void D_T_P2_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1;

  t1 = 6.0*xi;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 24.0;
  values[3] = -2.0+t1+2.0*eta;
  values[4] = -4.0+8.0*xi+6.0*eta;
  values[5] = -6.0+t1+12.0*eta;
}

// values of the derivatives in xi-xi direction
static void D_T_P2_2D_DeriveXiXi(double xi, double eta, 
                                       double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 12.0;
  values[4] = 6.0;
  values[5] = 2.0;
}

// values of the derivatives in xi-eta direction
static void D_T_P2_2D_DeriveXiEta(double xi, double eta, 
                                       double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 6.0;
  values[4] = 8.0;
  values[5] = 6.0;
}

// values of the derivatives in eta-eta direction
static void D_T_P2_2D_DeriveEtaEta(double xi, double eta, 
                                       double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 2.0;
  values[4] = 6.0;
  values[5] = 12.0;
}

// ***********************************************************************

TBaseFunct2D *BF_D_T_P2_2D_Obj = new TBaseFunct2D
        (6, BF_D_T_P2_2D, BFUnitTriangle, 
         D_T_P2_2D_Funct, D_T_P2_2D_DeriveXi,
         D_T_P2_2D_DeriveEta, D_T_P2_2D_DeriveXiXi,
         D_T_P2_2D_DeriveXiEta, D_T_P2_2D_DeriveEtaEta, 2, 2,
         0, NULL);
