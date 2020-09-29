#include <TimeConvDiff2D.h>
#include <MacroCell.h>

void ExampleFile()
{
  OutPut("Example: 2D_1D_Smmoth.h" << endl) ;

  TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=0;
  #define __SIMPATURS__

  #define __PBSConstT__
}

void InitialU1(double x, double y, double *values)
{

  values[0] = 0.;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}


// unit square (0,1)X(0,1)^2
void Exact(double x, double y,  double *values)
{
 double z, t;
//  const double Pi=3.14159265;
 
 z = TDatabase::ParamDB->REACTOR_P29;
 t = TDatabase::TimeDB->CURRENTTIME;
 double k = TDatabase::ParamDB->REACTOR_P28;
  
  values[0] = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[1] = Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[2] = -Pi*(exp(-k*t))*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[3] = -Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[4] = -3.*Pi*Pi*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z); ;  
}

void InitialCondition_Psd(double x, double y, double *values)
{
 double z, t;
//  double Pi =3.14159265;
 double k = TDatabase::ParamDB->REACTOR_P28;
 z = TDatabase::ParamDB->REACTOR_P29;
 t = 0;

  values[0] = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[1] = Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[2] = -Pi*(exp(-k*t))*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[3] = -Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[4] = -3.*Pi*Pi*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z); ;  
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
 double x, y, z, t;
//  double Pi =3.14159265;
 
 z = TDatabase::ParamDB->REACTOR_P29;
 t = TDatabase::TimeDB->CURRENTTIME;
 double k = TDatabase::ParamDB->REACTOR_P28;
 switch(BdComp)
 {
  case 0:
      x = Param;
      y =0;
  break;
  
  case 1:
      x = 1.;
      y = Param;            
  break;
  
  case 2:
      x = 1-Param;
      y =1.; 
  break;   
  
  case 3:
      x = 0.;
      y = 1.-Param;
  break;

  default: cout << "wrong boundary part number" << endl;
    break;
  }
    
  value = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
 }

void BoundCondition_LminLMax(BoundCond &cond_Lmin, BoundCond &cond_Lmax)
{
  cond_Lmin = DIRICHLET;
  cond_Lmax = DIRICHLET;
}

void BoundValue_LMin(double x, double y,  double *values)
{
 double z, t;
//  double Pi =3.14159265;
 
 z = TDatabase::ParamDB->REACTOR_P12;
 t = TDatabase::TimeDB->CURRENTTIME;
 double k = TDatabase::ParamDB->REACTOR_P28;
 
  values[0] = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[1] = Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[2] = -Pi*(exp(-k*t))*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[3] = -Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[4] = -3.*Pi*Pi*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z); ;  
}

void BoundValue_LMax(double x, double y,  double *values)
{
 double z, t;
//  double Pi =3.14159265;
 
 z = TDatabase::ParamDB->REACTOR_P13;
 t = TDatabase::TimeDB->CURRENTTIME;
 double k = TDatabase::ParamDB->REACTOR_P28;
 
  values[0] = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[1] = Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[2] = -Pi*(exp(-k*t))*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[3] = -Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[4] = -3.*Pi*Pi*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
}

// ========================================================================
// for scalar equations
// ========================================================================
void BilinearCoeffs_Psd(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps;
//   double Pi =3.14159265;
  int i;
  double *coeff, *param;
  double x, y, z, t;
  
  t = TDatabase::TimeDB->CURRENTTIME;
  double k = TDatabase::ParamDB->REACTOR_P28;
 
  if(TDatabase::ParamDB->REACTOR_P2)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P2;
  else
    eps = 0.;

  z = TDatabase::ParamDB->REACTOR_P29;
 
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
     }
    else
     {
      coeff[1] = 0.;  // u1
      coeff[2] = 0.;  // u2
     }
    coeff[3] = 0.;

    coeff[4] = (2.*eps*Pi*Pi - k)*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z); // f 
  }
}

// initial conditon
void InitialCondition_Psd_Intl(int N_Coord, double *X, double *values)
{
//  double Pi =3.14159265;
 double x, y, z, t = 0.;
  double k = TDatabase::ParamDB->REACTOR_P28;
  
   if(N_Coord!=3)
  {
   printf("N_Coord!=3 InitialCondition_psd_Intl %d!!!\n", N_Coord);  
#ifdef _MPI
     MPI_Finalize();
#endif  
     exit(0);   
  }
  x = X[0];
  y = X[1];    
  z = X[2];  
   
  values[0] = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[1] = Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[2] = -Pi*(exp(-k*t))*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[3] = -Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[4] = -3.*Pi*Pi*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);  
  
//  cout<< "Initial " << values[0] << endl;
  
  
}

void Exact_Psd_Intl(int N_Coord, double *X, double *values)
{

//  double Pi =3.14159265;
 double x, y, z;
 double t = TDatabase::TimeDB->CURRENTTIME;
 double k = TDatabase::ParamDB->REACTOR_P28;
  
 if(N_Coord!=3)
  {
   printf("N_Coord!=3 InitialCondition_psd_Intl %d!!!\n", N_Coord);  
#ifdef _MPI
     MPI_Finalize();
#endif  
     exit(0);   
  }
  x = X[0];
  y = X[1];    
  z = X[2];  
 
  values[0] = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[1] = Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[2] = -Pi*(exp(-k*t))*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[3] = -Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[4] = -3.*Pi*Pi*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z); 
}

// void BilinearCoeffs_Psd_Intl(int n_points, int N_Dim, double **Coords,
//                              double **parameters, double **coeffs)

void BilinearCoeffs_Psd_Intl(int n_points, int N_Dim, double **Coords,
                             double **parameters, double **coeffs)
{
  int i;
  double eps_L, *coeff;
  double x, y, z;
//   double Pi =3.14159265;
  double t = TDatabase::TimeDB->CURRENTTIME;
  double k = TDatabase::ParamDB->REACTOR_P28;
 
  if(TDatabase::ParamDB->REACTOR_P3)
    eps_L = 1.0/TDatabase::ParamDB->REACTOR_P3;
  else
    eps_L = 0.;

  
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    // param = parameters[i];

//     x = X[i];
//     y = Y[i];
//     z = Z[i];
    x = Coords[0][i];
    y = Coords[1][i];
    z = Coords[2][i];
    
    // diffusion
    coeff[0] = eps_L;
    
    // convection in z direction
    coeff[1] = 1;
    
    // reaction term, has to be negative for PBS coupled with C
    coeff[2] = 0;
    
    //rhs  
    coeff[3] =  (eps_L*Pi*Pi )*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z) 
                 -Pi* coeff[1]*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*sin(Pi*z);
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
    
//   for(i=0; i<N_V; i++)
//    cout<< i << " X[i] " << X[i] <<endl;
// 
// exit(0);

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


