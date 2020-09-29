// ***********************************************************************
// P1 element, discontinous, 2D, quadrilateral
// ***********************************************************************

// base function values
static void D_Q_P1_2D_Funct(double xi, double eta, double *values)
{
  values[0] = 1;
  values[1] = 3*xi;
  values[2] = 3*eta;
}

// values of the derivatives in xi direction
static void D_Q_P1_2D_DeriveXi(double xi, double eta, double *values)
{
  values[0] = 0;
  values[1] = 3;
  values[2] = 0;
}

// values of the derivatives in eta direction
static void D_Q_P1_2D_DeriveEta(double xi, double eta, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 3;
}

// values of the derivatives in xi-xi direction
static void D_Q_P1_2D_DeriveXiXi(double xi, double eta, 
                                       double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
}

// values of the derivatives in eta-eta direction
static void D_Q_P1_2D_DeriveEtaEta(double xi, double eta, 
                                       double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
}

// values of the derivatives in xi-eta direction
static void D_Q_P1_2D_DeriveXiEta(double xi, double eta, 
                                       double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
}
// ***********************************************************************

TBaseFunct2D *BF_D_Q_P1_2D_Obj = new TBaseFunct2D
        (3, BF_D_Q_P1_2D, BFUnitSquare, 
         D_Q_P1_2D_Funct, D_Q_P1_2D_DeriveXi,
         D_Q_P1_2D_DeriveEta, D_Q_P1_2D_DeriveXiXi,
         D_Q_P1_2D_DeriveXiEta, D_Q_P1_2D_DeriveEtaEta, 1, 1,
         0, NULL);
