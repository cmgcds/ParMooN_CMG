// ***********************************************************************
// P0 Raviart-Thomas vector element, nonconforming , 2D
// History:  02.09.2010 implementation (Alfonso)
// ***********************************************************************

// base function values
// vector function, orthoonal to edges, 
// function i has 
//  * flux 1 through edge i 
//  * flux 0 through other edges
static void N_T_RT0_2D_Funct(double xi, double eta, double *values)
{
  // first component
  values[0]= xi;
  values[1]= xi;
  values[2]= xi-1.;

   // second component
  values[3]= eta-1.;
  values[4]= eta;
  values[5]= eta;
}

// values of the derivatives in xi direction
static void N_T_RT0_2D_DeriveXi(double xi, double eta, double *values)
{
  // first component
  values[0]= 1.;
  values[1]= 1.;
  values[2]= 1.;

   // second component
  values[3]= 0.;
  values[4]= 0.;
  values[5]= 0.;
}

// values of the derivatives in eta direction
static void N_T_RT0_2D_DeriveEta(double xi, double eta, double *values)
{
  // first component
  values[0]= 0.;
  values[1]= 0.;
  values[2]= 0.;

   // second component
  values[3]= 1.;
  values[4]= 1.;
  values[5]= 1.;
}

// values of derivatives in xi-xi direction
static void N_T_RT0_2D_DeriveXiXi(double xi, double eta, double *values)
{
  // first component
  values[0]= 0.;
  values[1]= 0.;
  values[2]= 0.;
  
   // second component
  values[3]= 0.;
  values[4]= 0.;
  values[5]= 0;
}

// values of derivatives in eta-eta direction
static void N_T_RT0_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  // first component
  values[0]= 0;
  values[1]= 0;
  values[2]= 0;
  
   // second component
  values[3]= 0;
  values[4]= 0;
  values[5]= 0;
}

// values of derivatives in xi-eta direction
static void N_T_RT0_2D_DeriveXiEta(double xi, double eta, double *values)
{
  // first component
  values[0]= 0;
  values[1]= 0;
  values[2]= 0;
  
   // second component
  values[3]= 0;
  values[4]= 0;
  values[5]= 0;
}

// ***********************************************************************

TBaseFunct2D *BF_N_T_RT0_2D_Obj = new TBaseFunct2D
        (3, BF_N_T_RT0_2D, BFUnitTriangle, 
         N_T_RT0_2D_Funct, N_T_RT0_2D_DeriveXi,
         N_T_RT0_2D_DeriveEta, N_T_RT0_2D_DeriveXiXi,
         N_T_RT0_2D_DeriveXiEta, N_T_RT0_2D_DeriveEtaEta, 1, 1,
         0, NULL, 2);
