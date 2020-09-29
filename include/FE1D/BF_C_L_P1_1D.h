// ***********************************************************************
// P1 element, conforming, 1D
// ***********************************************************************

// base function values
static void C_L_P1_1D_Funct(double xi, double *values)
{
  values[0]=0.5*(1-xi);
  values[1]=0.5*(1+xi);
}

// values of the derivatives in xi direction
static void C_L_P1_1D_DeriveXi(double xi, double *values)
{
  values[0]=-0.5;
  values[1]= 0.5;
}

// values of the derivatives in xi-xi  direction
static void C_L_P1_1D_DeriveXiXi(double xi, double *values)
{
  values[0]=0;
  values[1]=0;
}

// ***********************************************************************

// TBaseFunct1D *BF_C_L_P1_1D_Obj = new TBaseFunct1D
//         (2, BF_C_L_P1_1D, C_L_P1_1D_Funct, C_L_P1_1D_DeriveXi, 
//          C_L_P1_1D_DeriveXiXi);
TBaseFunct1D *BF_C_L_P1_1D_Obj = new TBaseFunct1D
        (2, BF_C_L_P1_1D, C_L_P1_1D_Funct, C_L_P1_1D_DeriveXi, 
         C_L_P1_1D_DeriveXiXi, 1, 1);
