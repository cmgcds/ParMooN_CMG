// ***********************************************************************
// P2 element, conforming, 2D
// ***********************************************************************

// base function values
static void C_T_P2_2D_Funct(double xi, double eta, double *values)
{
  double t1 = xi+eta-1;

  values[0] = t1*(2*xi+2*eta-1);
  values[1] = -4*t1*xi;
  values[2] = xi*(2*xi-1);
  values[3] = -4*t1*eta;
  values[4] = 4*xi*eta;
  values[5] = eta*(2*eta-1);
}

// values of the derivatives in xi direction
static void C_T_P2_2D_DeriveXi(double xi, double eta, double *values)
{
  values[0] = 4*xi+4*eta-3;
  values[1] = -8*xi-4*eta+4;
  values[2] = 4*xi-1;
  values[3] = -4*eta;
  values[4] = 4*eta;
  values[5] = 0;
}

// values of the derivatives in eta direction
static void C_T_P2_2D_DeriveEta(double xi, double eta, double *values)
{
  values[0] = 4*xi+4*eta-3;
  values[1] = -4*xi;
  values[2] = 0;
  values[3] = -4*xi-8*eta+4;
  values[4] = 4*xi;
  values[5] = 4*eta-1;
}
// values of the derivatives in xi-xi  direction
static void C_T_P2_2D_DeriveXiXi(double xi, double eta, double *values)
{
  values[0] = 4;
  values[1] = -8;
  values[2] = 4;
  values[3] = 0;
  values[4] = 0;
  values[5] = 0;
}
// values of the derivatives in xi-eta direction
static void C_T_P2_2D_DeriveXiEta(double xi, double eta, double *values)
{
  values[0] = 4;
  values[1] = -4;
  values[2] = 0;
  values[3] = -4;
  values[4] = 4;
  values[5] = 0;
}
// values of the derivatives in eta-eta direction
static void C_T_P2_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  values[0] = 4;
  values[1] = 0;
  values[2] = 0;
  values[3] = -8;
  values[4] = 0;
  values[5] = 4;  
}

// ***********************************************************************

TBaseFunct2D *BF_C_T_P2_2D_Obj = new TBaseFunct2D
        (6, BF_C_T_P2_2D, BFUnitTriangle, 
         C_T_P2_2D_Funct, C_T_P2_2D_DeriveXi,
         C_T_P2_2D_DeriveEta, C_T_P2_2D_DeriveXiXi,
         C_T_P2_2D_DeriveXiEta, C_T_P2_2D_DeriveEtaEta, 2, 2,
         0, NULL);
