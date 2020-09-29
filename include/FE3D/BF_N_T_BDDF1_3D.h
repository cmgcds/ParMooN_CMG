// ***********************************************************************
// Brezzi-Douglas-Duran-Fortin element of first order on tetrahedra, 3D
// ***********************************************************************

static double N_T_BDDF1_3D_CM[144] = {
    0,-0,0,0,0,-0,0,0,-0,-1.666666667,0.3333333333,0.3333333333,
    0,0,0,-1.666666667,0.3333333333,0.3333333333,0,0,0,0,0,-0,
    -1.666666667,0.3333333333,0.3333333333,0,0,0,0,0,0,0,0,-0,
    -0.3333333333,1.666666667,-0.3333333333,-0.3333333333,-0.3333333333,1.666666667,-0.3333333333,1.666666667,-0.3333333333,1.666666667,-0.3333333333,-0.3333333333,
    0,0,0,2,0,-2,0,0,0,0,0,-0,
    2,-2,0,0,0,0,0,0,0,0,0,-0,
    -0,-0,0,0,-0,-0,0,-0,-0,2,-2,-0,
    -0.3333333333,-0.3333333333,1.666666667,1.666666667,-0.3333333333,-0.3333333333,1.666666667,-0.3333333333,-0.3333333333,-0.3333333333,1.666666667,-0.3333333333,
    2,0,-2,0,0,0,0,0,0,0,0,-0,
    0,0,0,0,0,0,0,0,0,2,0,-2,
    0,0,0,2,-2,0,0,0,0,0,0,-0,
    1.666666667,-0.3333333333,-0.3333333333,-0.3333333333,1.666666667,-0.3333333333,-0.3333333333,-0.3333333333,1.666666667,-0.3333333333,-0.3333333333,1.666666667
};

static void N_T_BDDF1_3D_Funct(double xi, double eta, double zeta,
                                                         double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={1,0,0,xi,0,0,eta,0,0,zeta,0,0};
  double mon_y[]={0,1,0,0,xi,0,0,eta,0,0,zeta,0};
  double mon_z[]={0,0,1,0,0,xi,0,0,eta,0,0,zeta};
  
  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_BDDF1_3D_DeriveXi(double xi, double eta, double zeta,
                                double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,1,0,0,0,0,0,0,0,0};
  double mon_y[]={0,0,0,0,1,0,0,0,0,0,0,0};
  double mon_z[]={0,0,0,0,0,1,0,0,0,0,0,0};
  
  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_BDDF1_3D_DeriveEta(double xi, double eta, double zeta,
                                 double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,1,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,1,0,0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,1,0,0,0};
  
  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_BDDF1_3D_DeriveZeta(double xi, double eta, double zeta,
                                  double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,1,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,1,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,1};
  
  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_BDDF1_3D_DeriveXiXi(double xi, double eta, double zeta,
                                  double *values)
{
  int nBF = 12; // number of basis functions
  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
}

static void N_T_BDDF1_3D_DeriveXiEta(double xi, double eta, double zeta,
                                   double *values)
{
  int nBF = 12; // number of basis functions
  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
}

static void N_T_BDDF1_3D_DeriveXiZeta(double xi, double eta, double zeta,
                                    double *values)
{
  int nBF = 12; // number of basis functions
  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
}

static void N_T_BDDF1_3D_DeriveEtaEta(double xi, double eta, double zeta,
                                    double *values)
{
  int nBF = 12; // number of basis functions
  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
}

static void N_T_BDDF1_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                                     double *values)
{
  int nBF = 12; // number of basis functions
  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
}

static void N_T_BDDF1_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                                      double *values)
{
  int nBF = 12; // number of basis functions
  memset(values, 0.0, 3*nBF*SizeOfDouble); // 3 is the space dimension
}

TBaseFunct3D *BF_N_T_BDDF1_3D_Obj =
new TBaseFunct3D(12, BF_N_T_BDDF1_3D, BFUnitTetrahedron,
                 N_T_BDDF1_3D_Funct, N_T_BDDF1_3D_DeriveXi,
                 N_T_BDDF1_3D_DeriveEta, N_T_BDDF1_3D_DeriveZeta,
                 N_T_BDDF1_3D_DeriveXiXi, N_T_BDDF1_3D_DeriveXiEta,
                 N_T_BDDF1_3D_DeriveXiZeta, N_T_BDDF1_3D_DeriveEtaEta,
                 N_T_BDDF1_3D_DeriveEtaZeta, N_T_BDDF1_3D_DeriveZetaZeta,
                 1, 1,
                 0, NULL, 3);
