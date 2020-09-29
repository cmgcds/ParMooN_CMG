#include <TimeConvDiff2D.h>
#include <MacroCell.h>

void ExampleFile()
{
  OutPut("Example: 2D_1D_ADI_UnitSrq.h" << endl) ;

  TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=1;
  #define __SIMPATURS__
}

void InitialU1(double x, double y, double *values)
{

  values[0] = 0.;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}


void InitialU2(double x, double y, double *values)
{
  values[0] = 0.;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}


// exact solution
void Exact(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}


void InitialCondition_Heat(double x, double y, double *values)
{
  double val;
  double T_D = TDatabase::ParamDB->REACTOR_P23;
  double T_W = TDatabase::ParamDB->REACTOR_P24;

  val = T_W/T_D;

  values[0] = val;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}


// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition_Heat(int i, double t, BoundCond &cond)
{

  switch(i)
  {
   case 0:
     cond = DIRICHLET;
   break;

   case 1:
     cond = NEUMANN;
   break;

   case 2:
   case 3:
   case 4:
   case 5:
     cond = DIRICHLET;
   break;

   default: 
      cout << "wrong boundary part number" << endl;
      exit(0);
   break;
  }

}

void BoundValue_Heat(int BdComp, double Param, double &value)
{
  double val;
  double T_D = TDatabase::ParamDB->REACTOR_P23;
  double T_W = TDatabase::ParamDB->REACTOR_P24;

  val = T_W/T_D;

  switch(BdComp)
  {
  case 1:
     value=0.;
  break;

  case 0:
  case 2:
  case 3:
  case 5:
    value=val;
  break;

  case 4:
      value = 1.;
  break;

  default:
     cout << "wrong boundary part number" << endl;
    break;
  }
}

// ========================================================================
// for scalar equations
// ========================================================================
void BilinearCoeffs_Heat(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps;
  int i;
  double *coeff, *param;
  double x, y;

  if(TDatabase::ParamDB->REACTOR_P0)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P0;
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
//       cout<< "coeff[0] eps " << eps << " u1 " << param[0]  << "u2 " << param[1] << endl;
     }
    else
     {
    coeff[1] = x;  // u1
    coeff[2] = -y;  // u2
     }
    coeff[3] = 0.;

    coeff[4] = 0.; // f
  }
}



void InitialCondition_Conc(double x, double y, double *values)
{

  values[0] = 1.;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}


// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition_Conc(int i, double t, BoundCond &cond)
{
  switch(i)
  {
   case 0:
   case 1:
   case 2:
   case 3:
   case 5:
     cond = NEUMANN;
   break;

   case 4:
     cond = DIRICHLET;
   break;

   default: 
      cout << "wrong boundary part number" << endl;
      exit(0);
   break;
  }
}

void BoundValue_Conc(int BdComp, double Param, double &value)
{

 switch(BdComp)
  {
  case 0:
  case 1:
  case 2:
  case 3:
  case 5:
     value=0.;
  break;

  case 4:
     value = 1.;
  break;

  default: cout << "wrong boundary part number" << endl;
    break;
  }
}


// ========================================================================
// for scalar equations
// ========================================================================
void BilinearCoeffs_Conc(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps;
//   double a=0.2, b=0., c=0.;
  int i;
  double *coeff, *param;
  double x, y;

  if(TDatabase::ParamDB->REACTOR_P1)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P1;
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
     }
    else
     {
      coeff[1] = x;  // u1
      coeff[2] = -y;  // u2
     }
    coeff[3] = 0.;

    coeff[4] = 0.; // f
  }
}

void InitialCondition_Psd(double x, double y, double *values)
{
  values[0] = 0.;
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
   case 5:
     cond = NEUMANN;
   break;

   case 4:
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
  case 5:
    value = 0.;
  break;
  case 4:
    value = 0.;
  break;

  default: cout << "wrong boundary part number" << endl;
    break;
  }
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
    coeff[1] = x;  // u1
    coeff[2] = -y;  // u2
     }
    coeff[3] = 0.;

    coeff[4] = 0.; // f
  }
}


// // initial conditon
// void PsdAtLmin(double x, double y, double z, double &value)
// {
// //   //testing the implementation
// //   double t = TDatabase::TimeDB->CURRENTTIME;
// //   value = x*t/5.;
// 
// //   value = value:
// }


// initial conditon
void InitialCondition_Psd_Intl(double x, double y, double z, double *values)
{

  //testing the implementation
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = 0;
//    cout<< " test example "<<  values[0] <<endl;
}

void BilinearCoeffs_Psd_Intl(int n_points, double *X, double *Y, double *Z,
        double **parameters, double **coeffs)
{
  int i;
  double eps, *coeff;                                  // *param;
  double x, y, z, c, a[3], b[3], s[3], h;
  double t = TDatabase::TimeDB->CURRENTTIME;

  b[0] = 0.1;
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
    coeff[1] = b[0];
    // reaction term
    coeff[2] = c;
    coeff[3] = 0.;
  }
}

void GetExampleFileData(BoundCondFunct2D **BoundaryConditions, BoundValueFunct2D **BoundValues, 
                        DoubleFunct2D **InitiaValues, CoeffFct2D **BilinearCoeffs, 
                        int &N_PBEqns, int &N_IndepntScalarEqns, int *Disctypes)
{

//   TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=1;


   N_IndepntScalarEqns = 2;
   N_PBEqns = 1;
   #define __PBS__

   BilinearCoeffs[0] = BilinearCoeffs_Heat;
   BilinearCoeffs[1] = BilinearCoeffs_Conc;
   BilinearCoeffs[2] = BilinearCoeffs_Psd;

   BoundaryConditions[0] = BoundCondition_Heat;
   BoundaryConditions[1] = BoundCondition_Conc;
   BoundaryConditions[2] = BoundCondition_Psd;

   BoundValues[0] = BoundValue_Heat;
   BoundValues[1] = BoundValue_Conc;
   BoundValues[2] = BoundValue_Psd;

   InitiaValues[0] = InitialCondition_Heat;
   InitiaValues[1] = InitialCondition_Conc;
   InitiaValues[2] = InitialCondition_Psd;

   Disctypes[0] = GALERKIN;
   Disctypes[1] = GALERKIN;
//    Disctypes[2] = SDFEM;
   Disctypes[2] = GALERKIN;
}


void Generate1DUniformMesh(TDomain *Domain, double Start, double End, int N_Cells)
{
  int i, j, N_V;
  int *Lines;
  double len, h, x, y, *X;
  TVertex **Vetrex;
  TJoint *Joint;
  TBaseCell  **CellTree;

  N_V = N_Cells+1;
  X = new double[N_V];

  for(i=0; i<N_V; i++)
   X[i] = 1. + (1. - Start)*(tanh(2.*(double(i)/double(N_Cells) - 1.)))/tanh(2.);

  X[0] = Start;
  X[N_V-1] = End;

//   for(i=0; i<N_V; i++)
//    cout<< " X[i] " << X[i] <<endl;

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


