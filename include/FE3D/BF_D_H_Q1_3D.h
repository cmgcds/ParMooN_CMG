// ***********************************************************************
// Q1 element, discontinuous, 3D
// ***********************************************************************

static void D_H_Q1_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
  values[0] = 1;
  values[1] = xi;
  values[2] = eta;
  values[3] = zeta;
  values[4] = xi*eta;
  values[5] = xi*zeta;
  values[6] = eta*zeta;
  values[7] = xi*eta*zeta;
}

static void D_H_Q1_3D_DeriveXi(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 1;
  values[2] = 0;
  values[3] = 0;
  values[4] = eta;
  values[5] = zeta;
  values[6] = 0;
  values[7] = eta*zeta;
}

static void D_H_Q1_3D_DeriveEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 1;
  values[3] = 0;
  values[4] = xi;
  values[5] = 0;
  values[6] = zeta;
  values[7] = xi*zeta;
}

static void D_H_Q1_3D_DeriveZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 1;
  values[4] = 0;
  values[5] = xi;
  values[6] = eta;
  values[7] = xi*eta;
}

static void D_H_Q1_3D_DeriveXiXi(double xi, double eta, double zeta,
                             double *values)
{
  memset(values, 0.0, 8*SizeOfDouble);
}

static void D_H_Q1_3D_DeriveXiEta(double xi, double eta, double zeta,
                             double *values)
{
  memset(values, 0.0, 8*SizeOfDouble);
  values[4] = 1;
  values[7] = zeta;
}

static void D_H_Q1_3D_DeriveXiZeta(double xi, double eta, double zeta,
                             double *values)
{
  memset(values, 0.0, 8*SizeOfDouble);
  values[5] = 1;
  values[7] = eta;
}

static void D_H_Q1_3D_DeriveEtaEta(double xi, double eta, double zeta,
                             double *values)
{
  memset(values, 0.0, 8*SizeOfDouble);
}

static void D_H_Q1_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                             double *values)
{
  memset(values, 0.0, 8*SizeOfDouble);
  values[6] = 1;
  values[7] = xi;
}

static void D_H_Q1_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                             double *values)
{
  memset(values, 0.0, 8*SizeOfDouble);
}

TBaseFunct3D *BF_D_H_Q1_3D_Obj = 
new TBaseFunct3D(8, BF_D_H_Q1_3D, BFUnitHexahedron, 
                 D_H_Q1_3D_Funct, D_H_Q1_3D_DeriveXi,
                 D_H_Q1_3D_DeriveEta, D_H_Q1_3D_DeriveZeta,
                 D_H_Q1_3D_DeriveXiXi, D_H_Q1_3D_DeriveXiEta,
                 D_H_Q1_3D_DeriveXiZeta, D_H_Q1_3D_DeriveEtaEta,
                 D_H_Q1_3D_DeriveEtaZeta, D_H_Q1_3D_DeriveZetaZeta,
                 3, 1,
                 0, NULL);
