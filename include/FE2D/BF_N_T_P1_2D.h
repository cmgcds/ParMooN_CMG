// ***********************************************************************
// P1 element, nonconforming, 2D
// ***********************************************************************

// base function values
static void N_T_P1_2D_Funct(double xi, double eta, double *values)
{
  values[0]=     -2*eta+1;
  values[1]= 2*xi+2*eta-1;
  values[2]=-2*xi      +1;
}

// values of the derivatives in xi direction
static void N_T_P1_2D_DeriveXi(double xi, double eta, double *values)
{
  values[0]= 0;
  values[1]= 2;
  values[2]=-2;
}

// values of the derivatives in eta direction
static void N_T_P1_2D_DeriveEta(double xi, double eta, double *values)
{
  values[0]=-2;
  values[1]= 2;
  values[2]= 0;
}

// values of the derivatives in xi-xi direction
static void N_T_P1_2D_DeriveXiXi(double xi, double eta, 
                                       double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
}

// values of the derivatives in eta-eta direction
static void N_T_P1_2D_DeriveEtaEta(double xi, double eta, 
                                       double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
}

// values of the derivatives in xi-eta direction
static void N_T_P1_2D_DeriveXiEta(double xi, double eta, 
                                       double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
}
// ***********************************************************************

TBaseFunct2D *BF_N_T_P1_2D_Obj = new TBaseFunct2D
        (3, BF_N_T_P1_2D, BFUnitTriangle, 
         N_T_P1_2D_Funct, N_T_P1_2D_DeriveXi,
         N_T_P1_2D_DeriveEta, N_T_P1_2D_DeriveXiXi,
         N_T_P1_2D_DeriveXiEta, N_T_P1_2D_DeriveEtaEta, 1, 1,
         0, NULL);
