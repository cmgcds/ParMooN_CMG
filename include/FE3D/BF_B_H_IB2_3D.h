// ***********************************************************************
// internal bubble of degree 2 (in the sense of Q2)
// ***********************************************************************

static void B_H_IB2_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
  values[0] = (1.0-xi)*(1.0-eta)*(1.0-zeta)*(1.0+xi)*(1.0+eta)*(1.0+zeta);
}

static void B_H_IB2_3D_DeriveXi(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = -2.0*xi+2.0*zeta*zeta*xi+2.0*eta*eta*xi-2.0*eta*eta*zeta*zeta*xi; 
}

static void B_H_IB2_3D_DeriveEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = -2.0*eta+2.0*eta*zeta*zeta+2.0*xi*xi*eta-2.0*zeta*zeta*xi*xi*eta;
}

static void B_H_IB2_3D_DeriveZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = -2.0*zeta+2.0*eta*eta*zeta+2.0*zeta*xi*xi-2.0*xi*xi*eta*eta*zeta;
}

static void B_H_IB2_3D_DeriveXiXi(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = -2.0*(-1.0+eta*eta)*(-1.0+zeta*zeta);
}

static void B_H_IB2_3D_DeriveXiEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = -4.0*eta*zeta*zeta*xi+4.0*eta*xi; 
}

static void B_H_IB2_3D_DeriveXiZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = -4.0*eta*eta*zeta*xi+4.0*zeta*xi; 
}

static void B_H_IB2_3D_DeriveEtaEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = -2.0*(-1.0+xi*xi)*(-1.0+zeta*zeta);
}

static void B_H_IB2_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 4.0*eta*zeta-4.0*zeta*xi*xi*eta; 
}

static void B_H_IB2_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = -2.0*(-1.0+xi*xi)*(-1.0+eta*eta);
}

TBaseFunct3D *BF_B_H_IB2_3D_Obj = 
new TBaseFunct3D(1, BF_B_H_IB2_3D, BFUnitHexahedron, 
                 B_H_IB2_3D_Funct,B_H_IB2_3D_DeriveXi,
                 B_H_IB2_3D_DeriveEta,B_H_IB2_3D_DeriveZeta,
                 B_H_IB2_3D_DeriveXiXi,B_H_IB2_3D_DeriveXiEta,
                 B_H_IB2_3D_DeriveXiZeta,B_H_IB2_3D_DeriveEtaEta,
                 B_H_IB2_3D_DeriveEtaZeta,B_H_IB2_3D_DeriveZetaZeta,
                 2, 2,
                 0, NULL);
