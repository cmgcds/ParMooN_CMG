// ***********************************************************************
// Raviart-Thomas element of zero-th order on tetrahedra, 3D
// ***********************************************************************

static void N_T_RT0_3D_Funct(double xi, double eta, double zeta,
                             double *values)
{
  // first component
  values[0] = xi;
  values[1] = xi;
  values[2] = xi;
  values[3] = xi - 1.0;
  
  // second component
  values[4] = eta;
  values[5] = eta - 1.0;
  values[6] = eta;
  values[7] = eta;
  
  // third component
  values[8] = zeta - 1.0;
  values[9] = zeta;
  values[10] = zeta;
  values[11] = zeta;
}

static void N_T_RT0_3D_DeriveXi(double xi, double eta, double zeta,
                                double *values)
{
  // first component
  values[0] = 1.0;
  values[1] = 1.0;
  values[2] = 1.0;
  values[3] = 1.0;
  
  // second component
  values[4] = 0;
  values[5] = 0;
  values[6] = 0;
  values[7] = 0;
  
  // third component
  values[8] = 0;
  values[9] = 0;
  values[10] = 0;
  values[11] = 0;
}

static void N_T_RT0_3D_DeriveEta(double xi, double eta, double zeta,
                                 double *values)
{
  // first component
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  
  // second component
  values[4] = 1.0;
  values[5] = 1.0;
  values[6] = 1.0;
  values[7] = 1.0;
  
  // third component
  values[8] = 0;
  values[9] = 0;
  values[10] = 0;
  values[11] = 0;
}

static void N_T_RT0_3D_DeriveZeta(double xi, double eta, double zeta,
                                  double *values)
{
  // first component
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  
  // second component
  values[4] = 0;
  values[5] = 0;
  values[6] = 0;
  values[7] = 0;
  
  // third component
  values[8] = 1.0;
  values[9] = 1.0;
  values[10] = 1.0;
  values[11] = 1.0;
}

static void N_T_RT0_3D_DeriveXiXi(double xi, double eta, double zeta,
                                  double *values)
{
  for(int i = 0; i < 12; i++)
    values[i] = 0;
}

static void N_T_RT0_3D_DeriveXiEta(double xi, double eta, double zeta,
                                   double *values)
{
  for(int i = 0; i < 12; i++)
    values[i] = 0;
}

static void N_T_RT0_3D_DeriveXiZeta(double xi, double eta, double zeta,
                                    double *values)
{
  for(int i = 0; i < 12; i++)
    values[i] = 0;
}

static void N_T_RT0_3D_DeriveEtaEta(double xi, double eta, double zeta,
                                    double *values)
{
  for(int i = 0; i < 12; i++)
    values[i] = 0;
}

static void N_T_RT0_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                                     double *values)
{
  for(int i = 0; i < 12; i++)
    values[i] = 0;
}

static void N_T_RT0_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                                      double *values)
{
  for(int i = 0; i < 12; i++)
    values[i] = 0;
}

TBaseFunct3D *BF_N_T_RT0_3D_Obj = 
new TBaseFunct3D(4, BF_N_T_RT0_3D, BFUnitTetrahedron, 
                 N_T_RT0_3D_Funct, N_T_RT0_3D_DeriveXi,
                 N_T_RT0_3D_DeriveEta, N_T_RT0_3D_DeriveZeta,
                 N_T_RT0_3D_DeriveXiXi, N_T_RT0_3D_DeriveXiEta,
                 N_T_RT0_3D_DeriveXiZeta, N_T_RT0_3D_DeriveEtaEta,
                 N_T_RT0_3D_DeriveEtaZeta, N_T_RT0_3D_DeriveZetaZeta,
                 1, 1,
                 0, NULL, 3);
