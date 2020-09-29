#include <TimeConvDiff2D.h>
#include <MacroCell.h>

void ExampleFile()
{
  OutPut("Example: SimPaTurS.h" << endl) ;
  #define __SIMPATURS__
  #define __PBS__
}

void InitialU1(double x, double y, double *values)
{

  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}


void InitialU2(double x, double y, double *values)
{

  values[0] = 0;
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
//   int rank;
  double T_W = 0.9864475690;

//   MPI_Comm_rank(MPI_COMM_WORLD, &rank);/**//**/
//   T_W = double(rank);

  values[0] = T_W;
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
   case 2: 
     cond = DIRICHLET;
   break;

   case 1: 
     cond = NEUMANN;
   break;

   case 3: 
   case 5: 
     cond = DIRICHLET;
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

void BoundValue_Heat(int BdComp, double Param, double &value)
{
  double val;
  double T_W = 0.9864475690;
  double T_In = 1.0;


  switch(BdComp)
  {
  case 0: 
  case 2:
     value=T_W;
  break;
  case 3:
  case 5: 
    value=T_W;
  break;
  case 1: 
    value=0.;
  break;
  case 4:
      value=T_In;
   break;
  default: cout << "wrong boundary part number" << endl;
    break;
  }
}

// ========================================================================
// for scalar equations
// ========================================================================
void BilinearCoeffs_Heat(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1./TDatabase::ParamDB->RE_NR;
  double a=0.2, b=0., c=0.;
  int i;
  double *coeff, *param;
  double x, y;

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
    coeff[1] = 0.;  // u1
    coeff[2] = 0.;  // u2
     }
    coeff[3] = c;

    coeff[4] = 0.; // f
  }
}



void InitialCondition_Conc(double x, double y, double *values)
{
//   int rank;
//   double T_W = 0.9864475690;
  double T_W = 0.5;


  values[0] = T_W;
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
   case 2: 
     cond = DIRICHLET;
   break;

   case 1: 
     cond = NEUMANN;
   break;

   case 3: 
   case 5: 
     cond = DIRICHLET;
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
  double val;
//   double T_W = 0.9864475690;
//   double T_In = 1.0;

  double T_W = -1;
  double T_In = 1.0;
  switch(BdComp)
  {
  case 0: 
  case 2:
     value=T_W;
  break;
  case 3:
  case 5: 
    value=T_W;
  break;
  case 1: 
    value=0.;
  break;
  case 4:
      value=T_In;
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
  double eps=1./TDatabase::ParamDB->RE_NR;
  double a=0.2, b=0., c=0.;
  int i;
  double *coeff, *param;
  double x, y;

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
    coeff[1] = 0.;  // u1
    coeff[2] = 0.;  // u2
     }
    coeff[3] = c;

    coeff[4] = 0.; // f
  }
}

void InitialCondition_Psd(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition_Psd(int i, double t, BoundCond &cond)
{
 cond = NEUMANN;
}

void BoundValue_Psd(int BdComp, double Param, double &value)
{
  value=0;
}

// ========================================================================
// for scalar equations
// ========================================================================
void BilinearCoeffs_Psd(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1./TDatabase::ParamDB->RE_NR;
  double a=0.2, b=0., c=0.;
  int i;
  double *coeff, *param;
  double x, y;

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
    coeff[1] = 0.;  // u1
    coeff[2] = 0.;  // u2
     }
    coeff[3] = c;

    coeff[4] = 0.; // f
  }
}


// initial conditon
void InitialCondition_Psd_Intl(double x, double y, double z, double *values)
{
    values[0] = 0;
}

void BilinearCoeffs_Psd_Intl(int n_points, double *X, double *Y, double *Z,
        double **parameters, double **coeffs)
{
  double eps = 1.0/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;                                  // *param;
  double x, y, z, c, a[3], b[3], s[3], h;
  double t = TDatabase::TimeDB->CURRENTTIME;
  
  b[0] = 1;
  b[1] = -2;
  b[2] = 2;  
  c = 0;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    // param = parameters[i];

    x = X[i];
    y = Y[i];
    z = Z[i];

    // diffusion
    coeff[0] = eps;
    // convection in x direction
    coeff[1] = b[0];
    // convection in y direction
    coeff[2] = b[1];
    // convection in z direction
    coeff[3] = b[2];
    // reaction
    coeff[4] = c;
     // rhs
    coeff[5] =0;
    // rhs from previous time step
    coeff[6] = 0;
  }
}



void GetExampleFileData(BoundCondFunct2D **BoundaryConditions, BoundValueFunct2D **BoundValues, 
                        DoubleFunct2D **InitiaValues, CoeffFct2D **BilinearCoeffs, 
                        int &N_PBEqns, int &N_IndepntScalarEqns)
{

//   TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=1;


   N_IndepntScalarEqns = 2;
   N_PBEqns = 1;

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
}



void Generate1DUniformMesh(TDomain *Domain, double Start, double End, int N_Cells)
{
  int i, j;
  int *Lines;
  double len, h, x, y;
  TVertex **Vetrex;
  TJoint *Joint;
  TBaseCell  **CellTree;

  len = End-Start;
  h = len/double(N_Cells);

  Lines = new int[2*N_Cells];
  Vetrex = new TVertex*[N_Cells+1]; 

  x=Start;
  y=0.;

  for(i=0; i<N_Cells; i++)
   {
    Lines[2*i]=i;
    Lines[2*i+1]=i+1;
    Vetrex[i] = new TVertex(x, y);
    x =double(i+1)*h;
   }

  x=End;
  Vetrex[N_Cells] = new TVertex(x, y);

  CellTree = new TBaseCell*[N_Cells];


   for (i=0;i<N_Cells;i++)
   {
//     Vetrex[ i ]->GetCoords(x, y);
//     cout<< " x " << x<< " y " << y<<endl;
    CellTree[i] = new TMacroCell(TDatabase::RefDescDB[S_Line], 0);
    CellTree[i]->SetVertex(0, Vetrex[ Lines[ 2*i       ]]);
    CellTree[i]->SetVertex(1, Vetrex[ Lines[ 2*i + 1]]);

    ((TMacroCell *) CellTree[i])->SetSubGridID(0);
//     ((TBaseCell *) SurfCellTree[i])->SetBd_Part(Bd_Part[i] );
   }

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







