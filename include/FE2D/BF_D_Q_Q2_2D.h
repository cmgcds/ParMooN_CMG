// ***********************************************************************
// Q2 element, discontinous, 2D, quadrilateral
// ***********************************************************************

// base function values
static void D_Q_Q2_2D_Funct(double xi, double eta, double *values)
{
  double t1, t6, t9, t13;

  t1 = xi*xi;
  t6 = -1.0+3.0*t1;
  t9 = eta*eta;
  t13 = -1.0+3.0*t9;

  values[0] = 1.0;
  values[1] = xi;
  values[2] = -1.0/2.0+3.0/2.0*t1;
  values[3] = eta;
  values[4] = xi*eta;
  values[5] = t6*eta/2.0;
  values[6] = -1.0/2.0+3.0/2.0*t9;
  values[7] = xi*t13/2.0;
  values[8] = t6*t13/4.0;
}

// values of the derivatives in xi direction
static void D_Q_Q2_2D_DeriveXi(double xi, double eta, double *values)
{
  double t4;

  t4 = eta*eta;

  values[0] = 0.0;
  values[1] = 1.0;
  values[2] = 3.0*xi;
  values[3] = 0.0;
  values[4] = eta;
  values[5] = 3.0*xi*eta;
  values[6] = 0.0;
  values[7] = -1.0/2.0+3.0/2.0*t4;
  values[8] = 3.0/2.0*xi*(-1.0+3.0*t4);
}

// values of the derivatives in eta direction
static void D_Q_Q2_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1;

  t1 = xi*xi;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 1.0;
  values[4] = xi;
  values[5] = -1.0/2.0+3.0/2.0*t1;
  values[6] = 3.0*eta;
  values[7] = 3.0*xi*eta;
  values[8] = 3.0/2.0*(-1.0+3.0*t1)*eta;
}

// values of the derivatives in xi-xi direction
static void D_Q_Q2_2D_DeriveXiXi(double xi, double eta, 
                                       double *values)
{
  double t2;

  t2 = eta*eta;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 3.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 3.0*eta;
  values[6] = 0.0;
  values[7] = 0.0;
  values[8] = -3.0/2.0+9.0/2.0*t2;
}

// values of the derivatives in xi-eta direction
static void D_Q_Q2_2D_DeriveXiEta(double xi, double eta, 
                                       double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 1.0;
  values[5] = 3.0*xi;
  values[6] = 0.0;
  values[7] = 3.0*eta;
  values[8] = 9.0*xi*eta;
}

// values of the derivatives in eta-eta direction
static void D_Q_Q2_2D_DeriveEtaEta(double xi, double eta, 
                                       double *values)
{
  double t2;

  t2 = xi*xi;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = 3.0;
  values[7] = 3.0*xi;
  values[8] = -3.0/2.0+9.0/2.0*t2;
}

// ***********************************************************************

TBaseFunct2D *BF_D_Q_Q2_2D_Obj = new TBaseFunct2D
        (9, BF_D_Q_Q2_2D, BFUnitSquare, 
         D_Q_Q2_2D_Funct, D_Q_Q2_2D_DeriveXi,
         D_Q_Q2_2D_DeriveEta, D_Q_Q2_2D_DeriveXiXi,
         D_Q_Q2_2D_DeriveXiEta, D_Q_Q2_2D_DeriveEtaEta, 2, 2,
         0, NULL);
