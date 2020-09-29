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
 
  
#ifdef _MPI
 MPI_Comm Comm;
 int rank;

 Comm = TDatabase::ParamDB->Comm;
 MPI_Comm_rank(Comm, &rank);

 if(rank==TDatabase::ParamDB->Par_P0)
#endif
 {  
  OutPut("Example: PBS.h " << endl);
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
//  double t = TDatabase::TimeDB->CURRENTTIME;
//  double k = 0.1;
//   
//   values[0] = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
//   values[1] = Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);
//   values[2] = -Pi*(exp(-k*t))*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
//   values[3] = -Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
//   values[4] = -3.*Pi*Pi*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);

  double t, T_int,T_W ;
   
  T_int =1.;  
  T_W = 2.;
  t = TDatabase::TimeDB->CURRENTTIME;
  
  values[0] = T_W  + (64.*(T_int - T_W)/(Pi*Pi*Pi))*(exp(-3.*Pi*Pi*t/4.))*cos(Pi*x/2.)*cos(Pi*y/2.)*cos(Pi*z/2.);
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;

}

// initial conditon
void InitialExact(double x, double y, double z, double *values)
{
  
//  double t = 0;
//  double k = 0.1;
//   
//   values[0] = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
//   values[1] = Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);
//   values[2] = -Pi*(exp(-k*t))*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
//   values[3] = -Pi*(exp(-k*t))*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
//   values[4] = -3.*Pi*Pi*(exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);  

  values[0] = 1.;
  values[1] = 0.;
  values[2] =  0.;
  values[3] =  0.;
  values[4] =  0.; 
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
//void BoundValue_c_A(int BdComp, double Param, double &value)
void BoundValue_temp(double x, double y, double z, double &value)
{
//  double t= TDatabase::TimeDB->CURRENTTIME;
//  double k = 0.1;
//  
//  value = (exp(-k*t))*sin(Pi*x)*cos(Pi*y)*cos(Pi*z); 
 
 value = 2.;
  
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
//     coeff[5] = (3.*eps*Pi*Pi - 0.1)*(exp(-0.1*t))*sin(Pi*x[i])*cos(Pi*y[i])*cos(Pi*z[i]); // f
    
     coeff[5] = 0;   
    coeff[6] = 0;  
  }
}

// ========================================================================
// definitions for the concentration
// ========================================================================

// initial condition
void InitialCondition_conc(double x, double y, double z, double *values)
{
  double eps = 1e-8, val, temp;  
  
  values[0] = 0;
    
  // inflow
  if ((fabs(x)<eps)&& (y>=0.33333333) && (y<=0.66666667) &&
      (z>=0.33333333) && (z<=0.66666667))
   {
     BoundValue_temp(x,y,z,temp);
     val = TDatabase::ParamDB->UREA_rho_sat_1+TDatabase::ParamDB->UREA_rho_sat_2*(temp - 273.15); 

     values[0] = val/(TDatabase::ParamDB->UREA_c_infty*TDatabase::ParamDB->UREA_m_mol);      
    //printf(" inflow temp \n" );
   } 
}

// kind of boundary condition (for FE space needed)
void BoundCondition_conc(double x, double y, double z, BoundCond &cond)
{
  
  double eps = 1e-8;
    
    
   cond = NEUMANN;
     
       // inflow
    if ((fabs(x)<eps)&& (y>=0.33333333) && (y<=0.66666667) &&
      (z>=0.33333333) && (z<=0.66666667))
    {
      cond = DIRICHLET;   
//       printf(" inflow BoundCondition_conc \n" );
    } 
    
    
}

// value of boundary condition
void BoundValue_conc(double x, double y, double z, double &value)
{
  double eps = 1e-8, val, temp;

  value = 0;
  // inflow
  if (TDatabase::TimeDB->CURRENTTIME >= TDatabase::TimeDB->T1)
  {
    if ((fabs(x)<eps)&& (y>=0.33333333) && (y<=0.66666667) &&
      (z>=0.33333333) && (z<=0.66666667))
      {
       BoundValue_temp(x,y,z,temp);
       val = TDatabase::ParamDB->UREA_rho_sat_1+TDatabase::ParamDB->UREA_rho_sat_2*(temp - 273.15); 

       value = val/(TDatabase::ParamDB->UREA_c_infty*TDatabase::ParamDB->UREA_m_mol);
       
//        printf("BoundValue_conc %f,\n", value);
      }
  }
}


// ========================================================================
// BilinearCoeffs for concentration 
// ========================================================================
void BilinearCoeffs_Conc(int n_points, double *x, double *y, double *z,
               double **parameters, double **coeffs)
{
  double eps;
  int i;
  double *coeff, *param;

  if(TDatabase::ParamDB->REACTOR_P1)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P1;
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
    coeff[4] = 0.;
    coeff[5] = 0.; // f     // rhs
    coeff[6] = 0.;  
    
  }
}


// ========================================================================
// definitions for the PSD
// ========================================================================

// initial condition
void InitialCondition_psd(double x, double y, double z, double *values)
{
  values[0] = 0;
}

void InitialCondition_psd_Intl(int N_Coord, double *X, double *values)
{
  values[0] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition_psd(double x, double y, double z, BoundCond &cond)
{
  double eps = 1e-8;
  
   cond = NEUMANN;
   
    if ((fabs(x)<eps)&& (y>=0.33333333) && (y<=0.66666667) &&
      (z>=0.33333333) && (z<=0.66666667))
     { cond = DIRICHLET; }
   
}

// value of boundary condition
void BoundValue_psd(double x, double y, double z, double &value)
{
  double eps = 1e-8, val, temp;
  double L = TDatabase::ParamDB->REACTOR_P29;
    
  // inflow
   if ((fabs(x)<eps)&& (y>=0.33333333) && (y<=0.66666667) &&
      (z>=0.33333333) && (z<=0.66666667))  
    {  value = InletPSD(L)/TDatabase::ParamDB->UREA_f_infty;
       // printf(" %e Inlet Urea %e \n", L,  value);  
    }
   else
    { value = 0; }
    
    
}


void BoundCondition_LminLMax(BoundCond &cond_Lmin, BoundCond &cond_Lmax)
{
  cond_Lmin = DIRICHLET;
//   cond_Lmin = NEUMANN;
  cond_Lmax = NEUMANN;
}

void BoundValue_LMin(double x, double y, double z, double *values)
 {
//    double t;
  
   // no nucleation, so f = 0 at L=0
   
//  if((x==0) && y>(1./3) &&  y<(2./3))
//   {
//   if((x==0) && y>(1./3) &&  y<(1./3. + 1./24.))
//    {
//     t = -6. + 12.*(y - 1./3.)*24.*(1. + 1./3. + 1./24. - y);
//     t = 1./(1. + exp(-t));
// //              cout << y << " t " << t << endl;
//     values[0] = t;
//     values[1] = 0.;
//     values[2] = 0.;    
//    }
//   else if((x==0) && y>(2./3 - 1./24.)  &&  y<(2./3))
//    {
//      t = 6. + 12.*(2./3. - 1./24. - y)*24.*(1. + 2./3. -y);
//      t = 1./(1. + exp(-t));
// //          cout << y << " t " << t << endl;
//     values[0] = t;
//     values[1] = 0.;
//     values[2] = 0.;    
//    }  
//   else
//    {
//     values[0] = 1.;
//     values[1] = 0.;
//     values[2] = 0.;
//    }
//   }
//  else
  {
    values[0] = 0.;
    values[1] = 0.;
    values[2] = 0.;
  } 

 }


void BoundValue_LMax(double x, double y,  double z,  double *values)
 {
    values[0] = 0.;
    values[1] = 0.;
    values[2] = 0.; 
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

  if(TDatabase::ParamDB->REACTOR_P2)
   { eps = 1.0/TDatabase::ParamDB->REACTOR_P2; }
  else
   { eps = 0.;}
  
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
     
    coeff[4] = 0; // reaction
    coeff[5] = 0.; // f // rhs
    coeff[6] = 0;  
  }
}


void BilinearCoeffs_Psd_Intl(int n_points, int N_Dim, double **Coords,
                             double **parameters, double **coeffs)
{
  int i;
  double eps, *coeff; // *param;
  double x, y, z, L, c, a[3], b[3], s[3], h;
  double t = TDatabase::TimeDB->CURRENTTIME;

  b[0] = -1e8;// negative, so that C will be taken from the PBS growth term
  c = 0; 

  if(TDatabase::ParamDB->REACTOR_P3)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P3;
  else
    eps = 0.;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    // param = parameters[i];

//     x = Coords[0][i];
//     y = Coords[1][i];
//     z = Coords[2][i];
//     L = Coords[3][i];    

//    printf("BilinearCoeffs_Psd_Intl  %f\n", L );
   
    // diffusion
    coeff[0] = eps;
    // convection in z direction
    coeff[1] = b[0];
    // reaction term
    coeff[2] = c;
    coeff[3] = 0.;
  }
}

void GetExampleFileData(BoundCondFunct3D **BoundaryConditions, BoundValueFunct3D **BoundValues, 
                        DoubleFunct3D **InitiaValues, CoeffFct3D **BilinearCoeffs, 
                        int &N_PBEqns, int &N_IndepntScalarEqns, int *Disctypes)
{
  
   N_IndepntScalarEqns = 1;
   N_PBEqns = 0;
   
//    #define __PBS__

   BilinearCoeffs[0] = BilinearCoeffs_Heat;
//    BilinearCoeffs[1] = BilinearCoeffs_Conc;
//    BilinearCoeffs[2] = BilinearCoeffs_Psd;

   BoundaryConditions[0] = BoundCondition_temp;
//    BoundaryConditions[1] = BoundCondition_conc;
//    BoundaryConditions[2] = BoundCondition_psd;

   BoundValues[0] = BoundValue_temp;
//    BoundValues[1] = BoundValue_conc;
//    BoundValues[2] = BoundValue_psd;

//    InitiaValues[0] = InitialCondition_temp;
   InitiaValues[0] = InitialExact; 
//    InitiaValues[0] = Exact;    
//    InitiaValues[1] = InitialCondition_conc;
//    InitiaValues[2] = InitialCondition_psd;

   Disctypes[0] = GALERKIN;
//    Disctypes[1] = GALERKIN;
// //    Disctypes[2] = SDFEM;
//    Disctypes[2] = GALERKIN;
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



