#include <TimeConvDiff2D.h>
#include <MacroCell.h>

void ExampleFile()
{
  OutPut("Example: 2D_1D_ConstT_UnitSqr_IntlOnly.h" << endl) ;

  TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=0;
  #define __SIMPATURS__

  #define __PBSConstT__
}


void Exact(double x, double y,  double *values)
{
  cout<< "Use Exact_Psd_Intl, see Example file " <<endl;
  exit(0);

}

void InitialCondition_Psd(double x, double y, double *values)
{
  cout<< "Use InitialCondition_Psd_Intl, see Example file " <<endl;
  exit(0);

}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition_Psd(int i, double t, BoundCond &cond)
{
  switch(i)
  {
   case 0:
   case 1:
   case 2:
   case 3:
     cond = DIRICHLET;
   break;

   default: 
      cout << "wrong boundary part number" << endl;
      exit(0);
   break;
  }
}

void BoundValue_Psd(int BdComp, double Param, double &value)
{
 switch(BdComp)
 {
  case 0:
  case 1:
  case 2:
  case 3:
    value = 0.;
  break;

  default: cout << "wrong boundary part number" << endl;
    break;
  }
}

void BoundCondition_LminLMax(BoundCond &cond_Lmin, BoundCond &cond_Lmax)
{
  cond_Lmin = DIRICHLET;
//   cond_Lmax = DIRICHLET;
  cond_Lmax = NEUMANN;
}


void BoundValue_LMin(double x, double y,  double *values)
 {
  values[0] = 1.;//u
  values[1] = 0.;//u_x
  values[2] = 0.; //u_xx
 }


void BoundValue_LMax(double x, double y,  double *values)
 {
  values[0] = 0.;  //u
  values[1] = 0.; //u_x
  values[2] = 0.; //u_xx
 }
// ========================================================================
// for scalar equations
// ========================================================================
void BilinearCoeffs_Psd(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
   int i;
  double *coeff, *param;
  double x, y;


  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = 0.;
    coeff[1] = 0.;  // u1
    coeff[2] = 0.;  // u2
    coeff[3] = 0.;
    coeff[4] = 0.; // f
  }
}

// initial conditon
void InitialCondition_Psd_Intl(double x, double y, double z, double *values)
{
 double t = 0.0;

  values[0] = (1. - z)*exp(-t*z*z);

}

void Exact_Psd_Intl(double x, double y, double z, double *values)
{

 double t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = (1. - z)*exp(-t*z*z);
}

void BilinearCoeffs_Psd_Intl(int n_points, double *X, double *Y, double *Z,
        double **parameters, double **coeffs)
{
  int i;
  double eps, *coeff;                                  // *param;
  double x, y, z, c, a[3], b[3], s[3], h;
  double t = TDatabase::TimeDB->CURRENTTIME;

  b[0] = 0;
  c = 0;

  if(TDatabase::ParamDB->REACTOR_P3)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P3;
  else
    eps = 0.;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    x = X[i];
    y = Y[i];
    z = Z[i];

    // diffusion
    coeff[0] = eps;
//     coeff[1] = 0.; // convection in z direction
    coeff[2] = 0.; // reaction term
//     coeff[3] = exp(-t*z*z)*(- z*z*(1.-z) +2.*t*(2.*t*z*z*z - 2.*t*z*z - 3.*z  - z  +1. )   );

//   with b=1
    coeff[0] = 0.;
    coeff[1] = 1.;
    coeff[3] = exp(-t*z*z)*(- z*z*(1.-z) 
                + coeff[0]*( 2.*t*(2.*t*z*z*z - 2.*t*z*z - 3.*z  - z  +1. ) +2.*t*z*z )
                - coeff[1]*( 2.*t*z  + 1. ));
  }
}

void GetExampleFileData(BoundCondFunct2D **BoundaryConditions, BoundValueFunct2D **BoundValues, 
                        DoubleFunct2D **InitiaValues, CoeffFct2D **BilinearCoeffs, 
                        int &N_PBEqns, int &N_IndepntScalarEqns, int *Disctypes)
{

  TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=0;


   N_IndepntScalarEqns = 0;
   N_PBEqns = 1;
   #define __PBS__

   BilinearCoeffs[0] = BilinearCoeffs_Psd;
   BoundaryConditions[0] = BoundCondition_Psd;
   BoundValues[0] = BoundValue_Psd;
   InitiaValues[0] = InitialCondition_Psd;
   Disctypes[0] = GALERKIN;

}

void Generate1DMesh(TDomain *Domain, double Start, double End, int N_Cells)
{
  int i, j, N_V;
  int *Lines;
  double len, h, x, y, *X;
  double hmin, hmax;  
  TVertex **Vetrex;
  TJoint *Joint;
  TBaseCell  **CellTree;

  N_V = N_Cells+1;
  X = new double[N_V];

  h = (End-Start)/(double)N_Cells;

  if(Start!=0. || End!=1.)
   {
    cout << "Domain should be (0,1) for this example " << endl;
    exit(0);
   }

  X[0] = Start;

  for(i=1; i<N_V; i++)
   X[i] = X[i-1] + (double)h;

  X[N_V-1] = End;

    hmin = 1.e8; hmax = -1.e8; 
    for(i=0; i<N_V-1; i++)
     {
      len = sqrt ((X[i+1] - X[i])*(X[i+1] - X[i]));
      if(len< hmin) hmin = len;
      if(len> hmax) hmax = len;
     }
     OutPut("L h_min : " << hmin << " L h_max : " << hmax << endl);

  Lines = new int[2*N_Cells];
  Vetrex = new TVertex*[N_V]; 

  y=0.;

  for(i=0; i<N_Cells; i++)
   {
    Lines[2*i]=i;
    Lines[2*i+1]=i+1;
    Vetrex[i] = new TVertex(X[i], y);
   }

  Vetrex[N_Cells] = new TVertex(X[N_V-1], y);

  CellTree = new TBaseCell*[N_Cells];

   for (i=0;i<N_Cells;i++)
   {
//     Vetrex[ i ]->GetCoords(x, y);
//     cout<< " x " << x<< " y " << y<<endl;
    CellTree[i] = new TMacroCell(TDatabase::RefDescDB[S_Line], 0);
    CellTree[i]->SetVertex(0, Vetrex[ Lines[ 2*i       ]]);
    CellTree[i]->SetVertex(1, Vetrex[ Lines[ 2*i + 1]]);
    ((TMacroCell *) CellTree[i])->SetSubGridID(0);
   }
//     Vetrex[ i ]->GetCoords(x, y);
//     cout<< " x " << x<< " y " << y<<endl;
//     exit(0);

   Domain->SetTreeInfo(CellTree, N_Cells);

   TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
   TDatabase::IteratorDB[It_LE]->SetParam(Domain);
   TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
   TDatabase::IteratorDB[It_Between]->SetParam(Domain);
   TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);

   // start joint(vertex)
   Joint = new TJointEqN(CellTree[0]);
   CellTree[0]->SetJoint(0, Joint);


   for(i=1;i<N_Cells;i++)
    {
     Joint = new TJointEqN(CellTree[i-1], CellTree[i]);

     CellTree[i-1]->SetJoint(1, Joint);
     CellTree[i]->SetJoint(0, Joint);
   } // for(i=0;i<N_Cells;i++)

   // end joint(vertex)
   Joint = new TJointEqN(CellTree[N_Cells-1]);
   CellTree[N_Cells-1]->SetJoint(1, Joint);

  delete []  Lines;
}


