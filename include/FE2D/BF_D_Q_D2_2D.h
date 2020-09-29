// ***********************************************************************
// P2 element, discontinous, 2D, quadrilateral
// ***********************************************************************

// base function values
static void D_Q_D2_2D_Funct(double xi, double eta, double *values)
{
  double t1, t3, t5, t7;

  t1 = xi*xi;
  t3 = 3.0/2.0*t1-1.0/2.0;
  t5 = eta*eta;
  t7 = 3.0/2.0*t5-1.0/2.0;

  values[0] = 1.0;
  values[1] = xi;
  values[2] = eta;
  values[3] = t3;
  values[4] = xi*eta;
  values[5] = t7;
  values[6] = t3*eta;
  values[7] = t7*xi;
  values[8] = 5.0/2.0*t1*xi*eta-5.0/2.0*xi*t5*eta;
}

// values of the derivatives in xi direction
static void D_Q_D2_2D_DeriveXi(double xi, double eta, double *values)
{
  double t4, t7;

  t4 = eta*eta;
  t7 = xi*xi;

  values[0] = 0.0;
  values[1] = 1.0;
  values[2] = 0.0;
  values[3] = 3.0*xi;
  values[4] = eta;
  values[5] = 0.0;
  values[6] = 3.0*xi*eta;
  values[7] = 3.0/2.0*t4-1.0/2.0;
  values[8] = 15.0/2.0*t7*eta-5.0/2.0*t4*eta;
}

// values of the derivatives in eta direction
static void D_Q_D2_2D_DeriveEta(double xi, double eta, double *values)
{
  double t2, t9;

  t2 = xi*xi;
  t9 = eta*eta;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 1.0;
  values[3] = 0.0;
  values[4] = xi;
  values[5] = 3.0*eta;
  values[6] = 3.0/2.0*t2-1.0/2.0;
  values[7] = 3.0*xi*eta;
  values[8] = 5.0/2.0*t2*xi-15.0/2.0*xi*t9;
}

// values of the derivatives in xi-xi direction
static void D_Q_D2_2D_DeriveXiXi(double xi, double eta, 
                                       double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 3.0;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = 3.0*eta;
  values[7] = 0.0;
  values[8] = 15.0*xi*eta;
}

// values of the derivatives in xi-eta direction
static void D_Q_D2_2D_DeriveXiEta(double xi, double eta, 
                                       double *values)
{
  double t3, t4;

  t3 = xi*xi;
  t4 = eta*eta;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 1.0;
  values[5] = 0.0;
  values[6] = 3.0*xi;
  values[7] = 3.0*eta;
  values[8] = 15.0/2.0*t3-15.0/2.0*t4;
}

// values of the derivatives in eta-eta direction
static void D_Q_D2_2D_DeriveEtaEta(double xi, double eta, 
                                       double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 3.0;
  values[6] = 0.0;
  values[7] = 3.0*xi;
  values[8] = -15.0*xi*eta;
}

// ***********************************************************************

TBaseFunct2D *BF_D_Q_D2_2D_Obj = new TBaseFunct2D
        (9, BF_D_Q_D2_2D, BFUnitSquare, 
         D_Q_D2_2D_Funct, D_Q_D2_2D_DeriveXi,
         D_Q_D2_2D_DeriveEta, D_Q_D2_2D_DeriveXiXi,
         D_Q_D2_2D_DeriveXiEta, D_Q_D2_2D_DeriveEtaEta, 3, 2,
         0, NULL);
