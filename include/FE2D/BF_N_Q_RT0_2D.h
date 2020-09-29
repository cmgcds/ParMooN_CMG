// ***********************************************************************
// Q0 Raviart-Thomas vector element, nonconforming , 2D
// History:  10.06.2010 implementation (Sashi)
//           01.07.2010 modification (Alfonso)
// ***********************************************************************

// base function values
// vector function, orthoonal to edges, 
// function i has 
//  * flux 1 through edge i 
//  * flux 0 through other edges
static void N_Q_RT0_2D_Funct(double xi, double eta, double *values)
{
  // first component
  values[0]= 0.;
  values[1]= 0.25*(xi+1.);
  values[2]= 0.;
  values[3]= 0.25*(xi-1.);

  // second component
  values[4]= 0.25*(eta-1.);
  values[5]= 0.;
  values[6]= 0.25*(eta+1.);
  values[7]= 0.; 
}

// values of the derivatives in xi direction
static void N_Q_RT0_2D_DeriveXi(double xi, double eta, double *values)
{
  // first component
  values[0]= 0.;
  values[1]= .25;
  values[2]= 0.;
  values[3]= .25;

   // second component
  values[4]= 0;
  values[5]= 0.;
  values[6]= 0;
  values[7]= 0.; 
}

// values of the derivatives in eta direction
static void N_Q_RT0_2D_DeriveEta(double xi, double eta, double *values)
{
  // first component
  values[0]= 0.;
  values[1]= 0.;
  values[2]= 0.;
  values[3]= 0.;
  
   // second component
  values[4]= .25;
  values[5]= 0;
  values[6]= .25;
  values[7]= 0; 
}

// values of derivatives in xi-xi direction
static void N_Q_RT0_2D_DeriveXiXi(double xi, double eta, double *values)
{
  // all second derivatives vanish
  memset(values, 0.0, 2*4*SizeOfDouble); // 2 is the space dimension
}

// values of derivatives in eta-eta direction
static void N_Q_RT0_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  // all second derivatives vanish
  memset(values, 0.0, 2*4*SizeOfDouble); // 2 is the space dimension
}

// values of derivatives in xi-eta direction
static void N_Q_RT0_2D_DeriveXiEta(double xi, double eta, double *values)
{
  // all second derivatives vanish
  memset(values, 0.0, 2*4*SizeOfDouble); // 2 is the space dimension
}

// ***********************************************************************

TBaseFunct2D *BF_N_Q_RT0_2D_Obj = new TBaseFunct2D
        (4, BF_N_Q_RT0_2D, BFUnitSquare, 
         N_Q_RT0_2D_Funct, N_Q_RT0_2D_DeriveXi,
         N_Q_RT0_2D_DeriveEta, N_Q_RT0_2D_DeriveXiXi,
         N_Q_RT0_2D_DeriveXiEta, N_Q_RT0_2D_DeriveEtaEta, 2, 1,
         0, NULL, 2);
