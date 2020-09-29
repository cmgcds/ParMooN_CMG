#include <TimeConvDiff2D.h>
#include <MacroCell.h>

void ExampleFile()
{
  OutPut("Example: 2D_1D_ConstT_UnitSrq.h" << endl) ;

  TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=0;
  #define __SIMPATURS__

  #define __PBSConstT__
}

// Chang et al., Numerical Heat Transfer, Part B vol 19, pp. 69-84, 1991
// exact solution for constant wall temperature for parallelepiped
// void Exact(double x, double y, double z, double *values)
// {
//  int i, j, k, N;
//  double t, a, s, T, T_int, T_w, temp;
//  double m, n, l;
//  double Phi =3.14159265;
// 
//  N= 100;
//  t = 0.2;
// 
//   a = 0.;
//   T_int =1.;
//   T_w = 2.;
// 
//   x = 0.5;
//   y = 0;
//   z = 0;
//   s = T_w;
//   for(i=1;i<N;i++)
//    for(j=1;j<N;j++)
//     for(k=1;k<N;k++)
//      {
//       m = (double)i;
//       n = (double)j;
//       l = (double)k;
//       a =(64.*(T_int - T_w)/(Phi*Phi*Phi*(2.*m-1.)*(2.*n-1.)*(2.*l-1.)))*sin((2.*m-1.)*Phi/2.)*sin((2.*n-1.)*Phi/2.)*sin((2.*l-1.)*Phi/2.);
//       temp = ((2.*m-1.)*Phi/2.)*((2.*m-1.)*Phi/2.) + ((2.*n-1.)*Phi/2.)*((2.*n-1.)*Phi/2.) + ((2.*l-1.)*Phi/2.)*((2.*l-1.)*Phi/2.);
//       s += a*(exp(-temp*t))*cos((2.*m-1.)*Phi*x/2.)*cos((2.*n-1.)*Phi*y/2.)*cos((2.*l-1.)*Phi*z/2.) ;
//      }
// 
//   values[0] = s;
//   values[1] = 0;
//   values[2] = 0;
//   values[3] = 0;
//   values[4] = 0;
// }

void Exact(double x, double y,  double *values)
{
 double t,  T_int, z;
 double Phi =3.14159265;

 z = TDatabase::ParamDB->REACTOR_P29;
 t = TDatabase::TimeDB->CURRENTTIME;

  T_int =1.;

  values[0] = (64.*T_int/(Phi*Phi*Phi))*(exp(-3.*Phi*Phi*t/4.))*cos(Phi*x/2.)*cos(Phi*y/2.)*cos(Phi*z/2.);
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void InitialCondition_Psd(double x, double y, double *values)
{
  values[0] = 1.;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
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
  cond_Lmax = DIRICHLET;
}


void BoundValue_LMin(double x, double y,  double *values)
 {
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.; 
 }


void BoundValue_LMax(double x, double y,  double *values)
 {
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;  
 }
// ========================================================================
// for scalar equations
// ========================================================================
void BilinearCoeffs_Psd(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps;
  int i;
  double *coeff, *param;
  double x, y;

  if(TDatabase::ParamDB->REACTOR_P2)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P2;
  else
    eps = 0.;


  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    if(TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
     {
      coeff[1] = param[0];  // u1
      coeff[2] = param[1];  // u2
      //cout<< "coeff[0] eps " << eps << " u1 " << param[0]  << "u2 " << param[1] << endl;
     }
    else
     {
      coeff[1] = 0.;  // u1
      coeff[2] = 0.;  // u2
     }
    coeff[3] = 0.;

    coeff[4] = 0.; // f
  }
}

// initial conditon
void InitialCondition_Psd_Intl(double x, double y, double z, double *values)
{
 double Phi =3.14159265;
 double t = 0.0;
 double T_int =1.;


  values[0] = (64.*T_int/(Phi*Phi*Phi))*(exp(-3.*Phi*Phi*t/4.))*cos(Phi*x/2.)*cos(Phi*y/2.)*cos(Phi*z/2.);
//   values[1] =

}

void Exact_Psd_Intl(double x, double y, double z, double *values)
{

 double Phi =3.14159265;
 double t = TDatabase::TimeDB->CURRENTTIME;
 double T_int =1.;

  values[0] = (64.*T_int/(Phi*Phi*Phi))*(exp(-3.*Phi*Phi*t/4.))*cos(Phi*x/2.)*cos(Phi*y/2.)*cos(Phi*z/2.);
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
    // param = parameters[i];

    x = X[i];
    y = Y[i];
    z = Z[i];

    // diffusion
    coeff[0] = eps;
    // convection in z direction
    coeff[1] = 0.;
    // reaction term
    coeff[2] = 0.;
    coeff[3] = 0.;
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


  if(Start!=-1. || End!=1.)
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


