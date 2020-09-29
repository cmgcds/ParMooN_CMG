// ***********************************************************************
// P00 element, conforming, 3D
// ***********************************************************************

static void C_T_P00_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
  values[0] = 0;
}

static void C_T_P00_3D_DeriveXi(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
}

static void C_T_P00_3D_DeriveEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
}

static void C_T_P00_3D_DeriveZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
}

static void C_T_P00_3D_DeriveXiXi(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
}

static void C_T_P00_3D_DeriveXiEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
}

static void C_T_P00_3D_DeriveXiZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
}

static void C_T_P00_3D_DeriveEtaEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
}

static void C_T_P00_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
}

static void C_T_P00_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
}

TBaseFunct3D *BF_C_T_P00_3D_Obj = 
new TBaseFunct3D(1, BF_C_T_P00_3D, BFUnitTetrahedron, 
                 C_T_P00_3D_Funct, C_T_P00_3D_DeriveXi,
                 C_T_P00_3D_DeriveEta, C_T_P00_3D_DeriveZeta,
                 C_T_P00_3D_DeriveXiXi, C_T_P00_3D_DeriveXiEta,
                 C_T_P00_3D_DeriveXiZeta, C_T_P00_3D_DeriveEtaEta,
                 C_T_P00_3D_DeriveEtaZeta, C_T_P00_3D_DeriveZetaZeta,
                 1, 1,
                 0, NULL);
