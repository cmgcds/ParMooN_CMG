// ***********************************************************************
// P3 element, conforming, 1D
// ***********************************************************************

// base function values
static void C_L_P3_1D_Funct(double xi, double *values)
{
  values[0]=-0.0625*(3*xi+1)*(3*xi-1)*(xi-1);
  values[1]= 0.5625*(xi+1)*(3*xi-1)*(xi-1);
  values[2]=-0.5625*(xi+1)*(3*xi+1)*(xi-1);
  values[3]= 0.0625*(xi+1)*(3*xi+1)*(3*xi-1);
}

// values of the derivatives in xi direction
static void C_L_P3_1D_DeriveXi(double xi, double *values)
{
  values[0]=0.0625*(xi*(-27*xi+18)+1);
  values[1]=0.0625*(xi*(81*xi-18)-27);
  values[2]=0.0625*(xi*(-81*xi-18)+27);
  values[3]=0.0625*(xi*(27*xi+18)+1);
}

// values of the derivatives in xi-xi  direction
static void C_L_P3_1D_DeriveXiXi(double xi, double *values)
{
  values[0]=0.125*(-27*xi+9);
  values[1]=0.125*( 81*xi-9);
  values[2]=0.125*(-81*xi-9);
  values[3]=0.125*( 27*xi+9);
}

// ***********************************************************************
// TBaseFunct1D *BF_C_L_P3_1D_Obj = new TBaseFunct1D
//         (4, BF_C_L_P3_1D, C_L_P3_1D_Funct, C_L_P3_1D_DeriveXi, 
//          C_L_P3_1D_DeriveXiXi);
TBaseFunct1D *BF_C_L_P3_1D_Obj = new TBaseFunct1D
        (4, BF_C_L_P3_1D, C_L_P3_1D_Funct, C_L_P3_1D_DeriveXi, 
         C_L_P3_1D_DeriveXiXi, 3, 3);
