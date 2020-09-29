// ***********************************************************************
// P0 element, conforming, 1D
// ***********************************************************************

// base function values
static void C_L_P0_1D_Funct(double xi, double *values)
{
  values[0]=1;
}

// values of the derivatives in xi direction
static void C_L_P0_1D_DeriveXi(double xi, double *values)
{
  values[0]=0;
}

// values of the derivatives in xi-xi  direction
static void C_L_P0_1D_DeriveXiXi(double xi, double *values)
{
  values[0]=0;
}

// ***********************************************************************
// TBaseFunct1D *BF_C_L_P0_1D_Obj = new TBaseFunct1D
//         (1, BF_C_L_P0_1D, C_L_P0_1D_Funct, C_L_P0_1D_DeriveXi, 
//          C_L_P0_1D_DeriveXiXi);

TBaseFunct1D *BF_C_L_P0_1D_Obj = new TBaseFunct1D
        (1, BF_C_L_P0_1D, C_L_P0_1D_Funct, C_L_P0_1D_DeriveXi, 
         C_L_P0_1D_DeriveXiXi, 0, 0);
