// ***********************************************************************
// P1 element, nonconforming, 3D
// ***********************************************************************

static void N_T_P1_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
  values[0] = 1-3*zeta;
  values[1] = 1-3*eta;
  values[2] = 3*(xi+eta+zeta)-2;
  values[3] = 1-3*xi;
}

static void N_T_P1_3D_DeriveXi(double xi, double eta, double zeta,
                             double *values)
{
  values[0] =  0;
  values[1] =  0;
  values[2] =  3;
  values[3] = -3;
}

static void N_T_P1_3D_DeriveEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] =  0;
  values[1] = -3;
  values[2] =  3;
  values[3] =  0;
}

static void N_T_P1_3D_DeriveZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = -3;
  values[1] =  0;
  values[2] =  3;
  values[3] =  0;
}

static void N_T_P1_3D_DeriveXiXi(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void N_T_P1_3D_DeriveXiEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void N_T_P1_3D_DeriveXiZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void N_T_P1_3D_DeriveEtaEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void N_T_P1_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void N_T_P1_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

TBaseFunct3D *BF_N_T_P1_3D_Obj = 
new TBaseFunct3D(4, BF_N_T_P1_3D, BFUnitTetrahedron, 
                 N_T_P1_3D_Funct, N_T_P1_3D_DeriveXi,
                 N_T_P1_3D_DeriveEta, N_T_P1_3D_DeriveZeta,
                 N_T_P1_3D_DeriveXiXi, N_T_P1_3D_DeriveXiEta,
                 N_T_P1_3D_DeriveXiZeta, N_T_P1_3D_DeriveEtaEta,
                 N_T_P1_3D_DeriveEtaZeta, N_T_P1_3D_DeriveZetaZeta,
                 1, 1,
                 0, NULL);
