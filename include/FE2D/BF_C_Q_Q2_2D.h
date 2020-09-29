// ***********************************************************************
// Q2 element, conforming, 2D
// ***********************************************************************

// base function values
static void C_Q_Q2_2D_Funct(double xi, double eta, double *values)
{
  double xi0=0.5*xi*(xi-1);
  double xi1=1-xi*xi;
  double xi2=0.5*xi*(xi+1);
  double eta0=0.5*eta*(eta-1);
  double eta1=1-eta*eta;
  double eta2=0.5*eta*(eta+1);

  values[0]=xi0*eta0;
  values[1]=xi1*eta0;
  values[2]=xi2*eta0;
  values[3]=xi0*eta1;
  values[4]=xi1*eta1;
  values[5]=xi2*eta1;
  values[6]=xi0*eta2;
  values[7]=xi1*eta2;
  values[8]=xi2*eta2;
}

// values of the derivatives in xi direction
static void C_Q_Q2_2D_DeriveXi(double xi, double eta, double *values)
{
  double xi0=xi-0.5;
  double xi1=-2*xi;
  double xi2=xi+0.5;
  double eta0=0.5*eta*(eta-1);
  double eta1=1-eta*eta;
  double eta2=0.5*eta*(eta+1);

  values[0]=xi0*eta0;
  values[1]=xi1*eta0;
  values[2]=xi2*eta0;
  values[3]=xi0*eta1;
  values[4]=xi1*eta1;
  values[5]=xi2*eta1;
  values[6]=xi0*eta2;
  values[7]=xi1*eta2;
  values[8]=xi2*eta2;
}

// values of the derivatives in eta direction
static void C_Q_Q2_2D_DeriveEta(double xi, double eta, double *values)
{
  double xi0=0.5*xi*(xi-1);
  double xi1=1-xi*xi;
  double xi2=0.5*xi*(xi+1);
  double eta0=eta-0.5;
  double eta1=-2*eta;
  double eta2=eta+0.5;

  values[0]=xi0*eta0;
  values[1]=xi1*eta0;
  values[2]=xi2*eta0;
  values[3]=xi0*eta1;
  values[4]=xi1*eta1;
  values[5]=xi2*eta1;
  values[6]=xi0*eta2;
  values[7]=xi1*eta2;
  values[8]=xi2*eta2;
}
// values of the derivatives in xi-xi  direction
static void C_Q_Q2_2D_DeriveXiXi(double xi, double eta, double *values)
{
  double eta0=0.5*eta*(eta-1);
  double eta1=1-eta*eta;
  double eta2=0.5*eta*(eta+1);

  values[0]=eta0;
  values[1]=-2*eta0;
  values[2]=eta0;
  values[3]=eta1;
  values[4]=-2*eta1;
  values[5]=eta1;
  values[6]=eta2;
  values[7]=-2*eta2;
  values[8]=eta2;
}
// values of the derivatives in xi-eta direction
static void C_Q_Q2_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double xi0=xi-0.5;
  double xi1=-2*xi;
  double xi2=xi+0.5;
  double eta0=eta-0.5;
  double eta1=-2*eta;
  double eta2=eta+0.5;

  values[0]=xi0*eta0;
  values[1]=xi1*eta0;
  values[2]=xi2*eta0;
  values[3]=xi0*eta1;
  values[4]=xi1*eta1;
  values[5]=xi2*eta1;
  values[6]=xi0*eta2;
  values[7]=xi1*eta2;
  values[8]=xi2*eta2;
}

// values of the derivatives in eta-eta direction
static void C_Q_Q2_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  double xi0=0.5*xi*(xi-1);
  double xi1=1-xi*xi;
  double xi2=0.5*xi*(xi+1);

  values[0]=xi0;
  values[1]=xi1;
  values[2]=xi2;
  values[3]=-2*xi0;
  values[4]=-2*xi1;
  values[5]=-2*xi2;
  values[6]=xi0;
  values[7]=xi1;
  values[8]=xi2; 
}


// ***********************************************************************

TBaseFunct2D *BF_C_Q_Q2_2D_Obj = new TBaseFunct2D
        (9, BF_C_Q_Q2_2D, BFUnitSquare, 
         C_Q_Q2_2D_Funct, C_Q_Q2_2D_DeriveXi,
         C_Q_Q2_2D_DeriveEta, C_Q_Q2_2D_DeriveXiXi,
         C_Q_Q2_2D_DeriveXiEta, C_Q_Q2_2D_DeriveEtaEta, 2, 2,
         0, NULL);
