#define __UREA__
#define __SIMPATURS__
  
#include <Urea_3d4d.h>
#include <MacroCell.h>

void ExampleFile()
{
  
  #define __PBSConstT__  
  
  
  // for velocity switch
  TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=0;

  TDatabase::ParamDB->N_CELL_LAYERS = 3;
  TDatabase::ParamDB->DRIFT_Z = 1;
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = 1356;
  TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE = FEM_FCT_LIN - 10000;   
  
#ifdef _MPI
 MPI_Comm Comm;
 int rank;

 Comm = TDatabase::ParamDB->Comm;
 MPI_Comm_rank(Comm, &rank);

 if(rank==TDatabase::ParamDB->Par_P0)
#endif
 {  
  OutPut("Example: PBS_UnitCubeSmoothT.h " << endl);
    
  OutPut("UREA_REACTION_DISC: " << TDatabase::ParamDB->UREA_REACTION_DISC << endl);
  OutPut("UREA_PB_DISC: " << TDatabase::ParamDB->UREA_PB_DISC << endl);
  OutPut("UREA_PB_DISC_STAB: " << TDatabase::ParamDB->UREA_PB_DISC_STAB<<endl);
  OutPut("UREA_SOLD_PARAMETER_TYPE: "<< TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE <<endl);
  OutPut("UREA_MODEL: " << TDatabase::ParamDB->UREA_MODEL << endl);
  OutPut("UREA_CONC_TOL: " << TDatabase::ParamDB->UREA_CONC_TOL << endl);
  OutPut("UREA_CONC_MAXIT: " << TDatabase::ParamDB->UREA_CONC_MAXIT << endl);

  OutPut("UREA_l_infty: " << TDatabase::ParamDB->UREA_l_infty <<endl);
  OutPut("UREA_u_infty: " << TDatabase::ParamDB->UREA_u_infty <<endl);
  OutPut("UREA_c_infty: " << TDatabase::ParamDB->UREA_c_infty <<endl);
  OutPut("UREA_temp_infty: " << TDatabase::ParamDB->UREA_temp_infty <<endl);
  OutPut("UREA_f_infty: " << TDatabase::ParamDB->UREA_f_infty<<endl);
  OutPut("UREA_nu: " << TDatabase::ParamDB->UREA_nu<<endl); 
  OutPut("UREA_rho: " << TDatabase::ParamDB->UREA_rho<<endl);
  OutPut("UREA_c_p: " << TDatabase::ParamDB->UREA_c_p<<endl);
  OutPut("UREA_lambda: " << TDatabase::ParamDB->UREA_lambda<<endl); 
  OutPut("UREA_D_P_0: " << TDatabase::ParamDB->UREA_D_P_0<<endl); 
  OutPut("UREA_D_P_MAX: " << TDatabase::ParamDB->UREA_D_P_MAX <<endl);
  OutPut("UREA_k_v: " << TDatabase::ParamDB->UREA_k_v<<endl); 
  OutPut("UREA_m_mol: " << TDatabase::ParamDB->UREA_m_mol<<endl);
  OutPut("UREA_D_J: " << TDatabase::ParamDB->UREA_D_J<<endl);
  OutPut("UREA_rho_d: " << TDatabase::ParamDB->UREA_rho_d <<endl);
  OutPut("UREA_k_g: " << TDatabase::ParamDB->UREA_k_g<<endl);
  OutPut("UREA_g: " << TDatabase::ParamDB->UREA_g <<endl);
  OutPut("UREA_rho_sat_1: " << TDatabase::ParamDB->UREA_rho_sat_1 <<endl);
  OutPut("UREA_rho_sat_2: " << TDatabase::ParamDB->UREA_rho_sat_2<<endl); 
  OutPut("UREA_beta_nuc: " << TDatabase::ParamDB->UREA_beta_nuc<<endl); 
  OutPut("UREA_alfa_nuc: " << TDatabase::ParamDB->UREA_alfa_nuc<<endl); 
  OutPut("UREA_INFLOW_SCALE: " << TDatabase::ParamDB->UREA_INFLOW_SCALE <<endl);
  OutPut("UREA_inflow_time: " << TDatabase::ParamDB->UREA_inflow_time <<endl);
 }
  // set some parameters
  //TDatabase::ParamDB->GRID_TYPE = 3;
  //OutPut("GRID_TYPE set to " << TDatabase::ParamDB->GRID_TYPE << endl);
}

// ========================================================================
// definitions for the temperature
// ========================================================================

void Exact(double x, double y, double z, double *values)
{
 double t = TDatabase::TimeDB->CURRENTTIME;
 double k = 0.1;
  
  values[0] = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[1] = Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[2] = -Pi*(exp(-k*t))*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[3] = -Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[4] = -3.*Pi*Pi*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
}

// initial conditon
void InitialExact(double x, double y, double z, double *values)
{ 
 double t = 0;
 double k = 0.1;
  
  values[0] = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[1] = Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[2] = -Pi*(exp(-k*t))*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[3] = -Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[4] = -3.*Pi*Pi*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);  
}


// initial conditon
void InitialCondition_temp(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;    
}


// kind of boundary condition (for FE space needed)
void BoundCondition_NSE(double x, double y, double z, BoundCond &cond)
{
  double eps = 1e-8;

  cond = DIRICHLET;

  if (fabs(x-210)<eps)
  {
       // outflow 
       cond = NEUMANN;
       //OutPut("neum");
       TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
  }
}

// kind of boundary condition (for FE space needed)
void BoundCondition_temp(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue_temp(double x, double y, double z, double &value)
{
 double t= TDatabase::TimeDB->CURRENTTIME;
 double k = 0.1;
 
 value = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);   
//      cout<< "value  " << value  << endl;
}


// ========================================================================
// BilinearCoeffs for Heat 
// ========================================================================
void BilinearCoeffs_Heat(int n_points, double *x, double *y, double *z,
               double **parameters, double **coeffs)
{
  double eps;
  int i;
  double *coeff, *param;
  double t = TDatabase::TimeDB->CURRENTTIME;
 
  if(TDatabase::ParamDB->REACTOR_P0)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P0;
  else
    eps = 0.;
  
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];    

    coeff[0] = eps;
    
    if(TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
     {
      coeff[1] = param[0];  // u1
      coeff[2] = param[1];  // u2
      coeff[3] = param[2];  // u2      
//       cout<< "coeff[0] eps " << eps << " u1 " << param[0]  << "u2 " << param[1] << endl;
     }
    else
     {
      coeff[1] = 0;  // u1
      coeff[2] = 0;  // u2
      coeff[3] = 0;  // u2      
     }
     
    // reaction
    coeff[4] = 0;
     // rhs
    coeff[5] = (3.*eps*Pi*Pi - 0.1)*(exp(-0.1*t))*sin(Pi*x[i])*cos(Pi*y[i])*cos(Pi*z[i]); // f
    coeff[6] = 0; 
//     cout<< "coeff[5]  " << coeff[5]  << endl;

  }
}
// ========================================================================
// definitions for the PSD
// ========================================================================

void Exact_Psd_Intl(double x, double y, double z, double l, double *values)
{
 double k = 0.1;
 double t= TDatabase::TimeDB->CURRENTTIME; 
 
  values[0] = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z)*sin(Pi*l);   
  values[1] = Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z)*sin(Pi*l);    
  values[2] = -Pi*(exp(-k*t))*sin(Pi*x)*sin(Pi*y)*cos(Pi*z)*sin(Pi*l);
  values[3] = -Pi*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*sin(Pi*z)*sin(Pi*l); 
  values[4] =  Pi*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z)*sin(Pi*l);   
  values[5] = -4*Pi*Pi*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z)*sin(Pi*l);    
}

void Exact_psd(int N_Coord, double *X, double *values)
{
 double t = TDatabase::TimeDB->CURRENTTIME; 
 
 if(N_Coord!=4)
  {
   printf("N_Coord!=4 InitialCondition_psd_Intl !!!\n");  
#ifdef _MPI
     MPI_Finalize();
#endif  
     exit(0);   
  }
 
  values[0] =  (exp(-0.1*t))*sin(Pi*X[0])*cos(Pi*X[1])*cos(Pi*X[2])*sin(Pi*X[3]);
}



// initial condition
void InitialCondition_psd(double x, double y, double z, double *values)
{
 double t = 0;  
 double k = 0.1;  
 double l = TDatabase::ParamDB->REACTOR_P29;
  
 values[0] = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z)*sin(Pi*l);
}

void InitialCondition_psd_Intl(int N_Coord, double *X, double *values)
{
 double t = 0;
  
 if(N_Coord!=4)
  {
   printf("N_Coord!=4 InitialCondition_psd_Intl !!!\n");  
#ifdef _MPI
     MPI_Finalize();
#endif  
     exit(0);   
  }
 
  values[0] = (exp(-0.1*t))*sin(Pi*X[0])*cos(Pi*X[1])*cos(Pi*X[2])*sin(Pi*X[3]);
}



// kind of boundary condition (for FE space needed)
void BoundCondition_psd(double x, double y, double z, BoundCond &cond)
{
   cond = DIRICHLET;  
}

// value of boundary condition
void BoundValue_psd(double x, double y, double z, double &value)
{
  double k = 0.1, temp;
  double l = TDatabase::ParamDB->REACTOR_P29;
  double t= TDatabase::TimeDB->CURRENTTIME;    

  value = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z)*sin(Pi*l);  
}


void BoundCondition_LminLMax(BoundCond &cond_Lmin, BoundCond &cond_Lmax)
{
  cond_Lmin = DIRICHLET;
//   cond_Lmin = NEUMANN;
  cond_Lmax = DIRICHLET;
}

void BoundValue_LMin(double x, double y, double z, double *values)
 {
  double k = 0.1;
  double l = TDatabase::ParamDB->REACTOR_P12;
  double t= TDatabase::TimeDB->CURRENTTIME;
 
  values[0] = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z)*sin(Pi*l); 
 }


void BoundValue_LMax(double x, double y,  double z,  double *values)
 {
  double k = 0.1;
  double l = TDatabase::ParamDB->REACTOR_P13;
  double t= TDatabase::TimeDB->CURRENTTIME;
 
  values[0] = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z)*sin(Pi*l); 
 }


// ========================================================================
// BilinearCoeffs for PSD 
// ========================================================================
void BilinearCoeffs_Psd(int n_points, double *x, double *y, double *z,
               double **parameters, double **coeffs)
{
  double eps;
  int i;
  double *coeff, *param;
  double t = TDatabase::TimeDB->CURRENTTIME;
  double l = TDatabase::ParamDB->REACTOR_P29;
  
  if(TDatabase::ParamDB->REACTOR_P0)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P0;
  else
    eps = 0.;
  
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];    

    coeff[0] = eps;
    
    if(TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
     {
      coeff[1] = param[0];  // u1
      coeff[2] = param[1];  // u2
      coeff[3] = param[2];  // u2      
     }
    else
     {
      coeff[1] = 0;  // u1
      coeff[2] = 0;  // u2
      coeff[3] = 0;  // u2      
     }
     
    // reaction
    coeff[4] = 0;
     // rhs
    coeff[5] = (3.*eps*Pi*Pi - 0.1)*(exp(-0.1*t))*sin(Pi*x[i])*cos(Pi*y[i])*cos(Pi*z[i])*sin(Pi*l); // f
    coeff[6] = 0;  
  }
}


void BilinearCoeffs_Psd_Intl(int n_points, int N_Dim, double **Coords,
                             double **parameters, double **coeffs)
{
  int i;
  double eps, *coeff; // *param;
  double x, y, z, L, c, a[3], b, s[3], h;
  double t = TDatabase::TimeDB->CURRENTTIME;

//   b = -1e8;// negative, so that C will be taken from the PBS growth term
  b = 0;
  c = 0; 

  if(TDatabase::ParamDB->REACTOR_P3)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P3;
  else
    eps = 0.;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    x = Coords[0][i];
    y = Coords[1][i];
    z = Coords[2][i];
    L = Coords[3][i];    
   
    // diffusion
    coeff[0] = eps;   
    
    // convection in z direction
    coeff[1] = b;
    // reaction term
    coeff[2] = c;
    coeff[3] = eps*Pi*Pi*(exp(-0.1*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z)*sin(Pi*L); // rhs
  }
}

void GetExampleFileData(BoundCondFunct3D **BoundaryConditions, BoundValueFunct3D **BoundValues, 
                        DoubleFunct3D **InitiaValues, CoeffFct3D **BilinearCoeffs, 
                        int &N_PBEqns, int &N_IndepntScalarEqns, int *Disctypes)
{
  
   N_IndepntScalarEqns = 0;
   N_PBEqns = 1;
   
     #define __PBS__

//    BilinearCoeffs[0] = BilinearCoeffs_Heat;
   BilinearCoeffs[0] = BilinearCoeffs_Psd;

//    BoundaryConditions[0] = BoundCondition_temp;
   BoundaryConditions[0] = BoundCondition_psd;

//    BoundValues[0] = BoundValue_temp;
   BoundValues[0] = BoundValue_psd;

//    InitiaValues[0] = InitialExact; 
   InitiaValues[0] = InitialCondition_psd;

   Disctypes[0] = GALERKIN;
 
}

void Generate1DMesh(TDomain *Domain, double Start, double End, int N_Cells)
{
  int i, j, N_V;
  int *Lines;
  double len, h, x, y, z, *X;
  TVertex **Vetrex;
  TJoint *Joint;
  TBaseCell  **CellTree;

  N_V = N_Cells+1;
  X = new double[N_V];

  for(i=0; i<N_V; i++)
   X[i] = 1. + (1. - Start)*(tanh(2.75*(double(i)/double(N_Cells) - 1.)))/tanh(2.75);

//   h = (End - Start)/N_Cells;
//   for(i=1; i<N_V; i++)
//    X[i] =  h*(double)i;

  X[0] = Start;
  X[N_V-1] = End;

//   for(i=0; i<N_V; i++)
//    cout<< " X[i] " << X[i] <<endl;

  Lines = new int[2*N_Cells];
  Vetrex = new TVertex*[N_V]; 

  y=0.;
  z=0.;
  for(i=0; i<N_Cells; i++)
   {
    Lines[2*i]=i;
    Lines[2*i+1]=i+1;
    Vetrex[i] = new TVertex(X[i], y, z);
   }

  Vetrex[N_Cells] = new TVertex(X[N_V-1], y, z);

  CellTree = new TBaseCell*[N_Cells];

   for (i=0;i<N_Cells;i++)
   {
//     Vetrex[ i ]->GetCoords(x, y);
//     cout<< " x " << x<< " y " << y<<endl;
    CellTree[i] = new TMacroCell(TDatabase::RefDescDB[S_Line], 0);
    CellTree[i]->SetVertex(0, Vetrex[ Lines[ 2*i       ]]);
    CellTree[i]->SetVertex(1, Vetrex[ Lines[ 2*i + 1]]);
    ((TMacroCell *) CellTree[i])->SetSubGridID(0);

//  cout<< " x " <<CellTree[i]->GetN_Edges()<<endl;;
//  cout<< " x " <<TDatabase::RefDescDB[S_Line]->GetN_OrigEdges()<<endl;;

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



