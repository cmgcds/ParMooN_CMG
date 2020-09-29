// ***********************************************************************
// P1 element, discontinuous, 3D Tetrahedra
// 
// Author:     Sashikumaar Ganesan
//
// ***********************************************************************

static void D_T_P1_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
   values[0] =  6.0;
   values[1] = -1.0+2.0*xi+eta+zeta;
   values[2] = -1.0+xi+2.0*eta+zeta;
   values[3] = -1.0+xi+eta+2.0*zeta;
}

static void D_T_P1_3D_DeriveXi(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 2;
  values[2] = 1;
  values[3] = 1;
}

static void D_T_P1_3D_DeriveEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 1;
  values[2] = 2;
  values[3] = 1;
}

static void D_T_P1_3D_DeriveZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 1;
  values[2] = 1;
  values[3] = 2;
}

static void D_T_P1_3D_DeriveXiXi(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void D_T_P1_3D_DeriveXiEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void D_T_P1_3D_DeriveXiZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void D_T_P1_3D_DeriveEtaEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void D_T_P1_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void D_T_P1_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

TBaseFunct3D *BF_D_T_P1_3D_Obj =
new TBaseFunct3D(4, BF_D_T_P1_3D, BFUnitTetrahedron,
                 D_T_P1_3D_Funct, D_T_P1_3D_DeriveXi,
                 D_T_P1_3D_DeriveEta, D_T_P1_3D_DeriveZeta,
                 D_T_P1_3D_DeriveXiXi, D_T_P1_3D_DeriveXiEta,
                 D_T_P1_3D_DeriveXiZeta, D_T_P1_3D_DeriveEtaEta,
                 D_T_P1_3D_DeriveEtaZeta, D_T_P1_3D_DeriveZetaZeta,
                 1, 1,
		 0, NULL);
