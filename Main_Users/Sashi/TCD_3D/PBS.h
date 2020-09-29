#define __UREA__
#define __SIMPATURS__
  
#include <Urea_3d4d.h>
#include <MacroCell.h>

void ExampleFile()
{
  // for velocity switch
  TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=1;

  TDatabase::ParamDB->N_CELL_LAYERS = 3;
  TDatabase::ParamDB->DRIFT_Z = 1;
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = 1356;
 
  TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE = FEM_FCT_LIN - 10000;  

  TDatabase::ParamDB->UREA_temp_infty = TDatabase::ParamDB->REACTOR_P23;
  TDatabase::ParamDB->UREA_c_infty = TDatabase::ParamDB->REACTOR_P25;
  TDatabase::ParamDB->UREA_f_infty = TDatabase::ParamDB->REACTOR_P20;
  TDatabase::ParamDB->UREA_rho_sat_1 = 35.36400;
  TDatabase::ParamDB->UREA_rho_sat_2= 1.305000;   
  TDatabase::ParamDB->UREA_MODEL=2;
  
  TDatabase::ParamDB->UREA_AGGR_SPATIAL =  3;
  TDatabase::ParamDB->UREA_AGGR_BROWNIAN = 1;
  TDatabase::ParamDB->UREA_AGGR_POL_ORDER = 1;
  TDatabase::ParamDB->UREA_AGGR_BROWNIAN_TEMP = 1; 
  TDatabase::ParamDB->UREA_AGGR_SHEAR_FACTOR_TYPE = 1;
  TDatabase::ParamDB->UREA_AGGR_BROWNIAN_SCAL = 2.0e5;  // ???  
  TDatabase::ParamDB->UREA_AGGR_SHEAR_FACTOR = 0.01;   // ???
  
  TDatabase::ParamDB->UREA_D_P_MAX =  5.0e-3;  
//   TDatabase::ParamDB->UREA_m_mol = 1.;
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
    
  OutPut("UREA_AGGR_SPATIAL  " << TDatabase::ParamDB->UREA_AGGR_SPATIAL <<endl);
  OutPut("UREA_AGGR_BROWNIAN  " << TDatabase::ParamDB->UREA_AGGR_BROWNIAN  <<endl);
  OutPut("UREA_AGGR_POL_ORDER  " << TDatabase::ParamDB->UREA_AGGR_POL_ORDER <<endl);
  
  OutPut("UREA_AGGR_BROWNIAN_TEMP  " << TDatabase::ParamDB->UREA_AGGR_BROWNIAN_TEMP  <<endl);
  OutPut("UREA_AGGR_BROWNIAN_SCAL  " << TDatabase::ParamDB->UREA_AGGR_BROWNIAN_SCAL <<endl);  
  
  OutPut("UREA_AGGR_SHEAR_FACTOR_TYPE  " <<  TDatabase::ParamDB->UREA_AGGR_SHEAR_FACTOR_TYPE <<endl);
  OutPut("UREA_AGGR_SHEAR_FACTOR  " << TDatabase::ParamDB->UREA_AGGR_SHEAR_FACTOR <<endl);     
 }
  // set some parameters
  //TDatabase::ParamDB->GRID_TYPE = 3;
  //OutPut("GRID_TYPE set to " << TDatabase::ParamDB->GRID_TYPE << endl);
}


int Is_InletBD(double x, double y, double z)
{
  int value=0;
  double eps = 1e-8;
 
    if( (fabs(x)<eps)&& (y>=0.33333333) && (y<=0.66666667) &&
        (z>=0.33333333) && (z<=0.66666667) )
     { value = 1;}
    
  return value;
}



void Exact(double x, double y, double z, double *values)
{ 
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void Exact_Psd_Intl(double x, double y, double z, double l, double *values)
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

// ========================================================================
// definitions for the temperature
// ========================================================================

// initial conditon
void InitialCondition_temp(double x, double y, double z, double *values)
{    
  values[0] = TDatabase::ParamDB->REACTOR_P23/TDatabase::ParamDB->UREA_temp_infty;
}

// kind of boundary condition (for FE space needed)
void BoundCondition_temp(double x, double y, double z, BoundCond &cond)
{
  double eps = 1e-8;

  cond = DIRICHLET;

  if (fabs(x-210)<eps)  // outflow 
   {      
    cond = NEUMANN;
   }
}

// value of boundary condition
//void BoundValue_c_A(int BdComp, double Param, double &value)
void BoundValue_temp(double x, double y, double z, double &value)
{
  double eps = 1e-8;
  double T_D = TDatabase::ParamDB->REACTOR_P23;
  double T_W = TDatabase::ParamDB->REACTOR_P24;
  
  value = 0;

    if( Is_InletBD(x, y, z) )
     { value = T_D/TDatabase::ParamDB->UREA_temp_infty; } 
    else
     {
      // outlet
      if ((fabs(x-210)<eps)&&(fabs(y)>eps)&&(fabs(z)>eps)&&(fabs(y-1.0)>eps)&&(fabs(z-1.0)>eps))
      {  value = 0; }
      else // wall
      {  value = T_W/TDatabase::ParamDB->UREA_temp_infty; }
     }

}

// ========================================================================
// BilinearCoeffs for Heat 
// ========================================================================
void BilinearCoeffs_Heat(int n_points, double *x, double *y, double *z,
               double **parameters, double **coeffs)
{
  int i;
  double eps, *coeff, *param;

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
      coeff[3] = param[2];  // u3      
//       cout<< "coeff[0] eps " << eps << " u1 " << param[0]  << "u2 " << param[1] << endl;
     }
    else
     {
      coeff[1] = 0;  // u1
      coeff[2] = 0;  // u2
      coeff[3] = 0;  // u2      
     }
    
    coeff[4] = 0;
    coeff[5] = 0.;
    coeff[6] = 0;  
  }
}

// ========================================================================
// definitions for the concentration
// ========================================================================
// initial condition
void InitialCondition_conc(double x, double y, double z, double *values)
{
  double val;
  
  val  = TDatabase::ParamDB->REACTOR_P23;
  val = TDatabase::ParamDB->UREA_rho_sat_1 + TDatabase::ParamDB->UREA_rho_sat_2*(val - 273.15); 
  val /=TDatabase::ParamDB->UREA_m_mol;
//        printf("InitialCondition_conc %f,\n", val);    
  values[0] = val/TDatabase::ParamDB->UREA_c_infty;
}

// kind of boundary condition (for FE space needed)
void BoundCondition_conc(double x, double y, double z, BoundCond &cond)
{
    
   cond = NEUMANN;
     
    // inflow
    if( Is_InletBD(x, y, z) )
     { cond = DIRICHLET; }   
    
}

// value of boundary condition
void BoundValue_conc(double x, double y, double z, double &value)
{
  double eps = 1e-8, val, temp;

  value = 0;
  
  // inflow
  if (TDatabase::TimeDB->CURRENTTIME >= TDatabase::TimeDB->T1)
  {
     if( Is_InletBD(x, y, z) )
      {
       BoundValue_temp(x,y,z,temp);
       temp *=TDatabase::ParamDB->UREA_temp_infty;
       val = TDatabase::ParamDB->UREA_rho_sat_1 + TDatabase::ParamDB->UREA_rho_sat_2*(temp - 273.15); 
       val /=TDatabase::ParamDB->UREA_m_mol;
       
       value = val/TDatabase::ParamDB->UREA_c_infty;       
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
  int i;
  double eps, *coeff, *param;

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
      coeff[3] = param[2];  // u3      
//       cout<< "coeff[0] eps " << eps << " u1 " << param[0]  << "u2 " << param[1] << endl;
     }
    else
     {
      coeff[1] = 0;  // u1
      coeff[2] = 0;  // u2
      coeff[3] = 0;  // u2      
     }
     
    coeff[4] = 0;
    coeff[5] = 0.;
    coeff[6] = 0;  
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
   
   if( Is_InletBD(x, y, z) )
    { cond = DIRICHLET; }
   
}


double Inlet_PSD(double L)  
{
 int i;
 double value, m, t, eps = 1e-8, V_inj;
 double L_max =  TDatabase::ParamDB->UREA_D_P_MAX;
 
  double inlet_coord[100] ={2.500000e-06,1.356050e-05,3.068150e-05,4.780300e-05,6.492450e-05,8.204550e-05,9.916650e-05,1.162875e-04,1.334090e-04,1.505305e-04,1.676515e-04,1.847725e-04,2.018940e-04,2.190155e-04,2.361365e-04,2.532575e-04,2.703785e-04,2.875000e-04,3.046215e-04,3.217425e-04,3.388635e-04,3.559845e-04,3.731060e-04,3.902275e-04,4.073485e-04,4.244695e-04,4.415910e-04,4.587125e-04,4.758335e-04,4.929545e-04,5.100755e-04,5.271970e-04,5.443185e-04,5.614395e-04,5.785605e-04,5.956815e-04,6.128030e-04,6.299245e-04,6.470455e-04,6.641665e-04,6.812875e-04,6.984090e-04,7.155305e-04,7.326515e-04,7.497725e-04,7.668940e-04,7.840155e-04,8.011365e-04,8.182575e-04,8.353785e-04,8.525000e-04,8.696215e-04,8.867425e-04,9.038635e-04,9.209845e-04,9.381060e-04,9.552275e-04,9.723485e-04,9.894695e-04,1.006591e-03,1.023712e-03,1.040833e-03,1.057955e-03,1.075075e-03,1.092197e-03,1.109318e-03,1.126439e-03,1.143561e-03,1.160682e-03,1.177803e-03,1.194925e-03,1.212045e-03,1.229167e-03,1.246287e-03,1.263409e-03,1.280530e-03,1.297651e-03,1.314773e-03,1.331894e-03,1.349016e-03,1.366137e-03,1.383257e-03,1.400379e-03,1.417500e-03,1.434622e-03,1.451742e-03,1.468863e-03,1.485985e-03,1.503106e-03,1.520228e-03,1.537348e-03,1.554470e-03,1.571591e-03,1.588713e-03,1.605834e-03,1.622954e-03,1.640075e-03,1.657197e-03,1.674318e-03,1.691440e-03};

  double inlet_f_L_seed[100] = {0.000000e+00,1.539144e+09,2.082543e+09,1.056886e+09,1.032820e+09,9.174217e+08,7.339395e+08,6.383108e+08,5.630089e+08,4.776040e+08,3.952118e+08,3.244646e+08,2.668657e+08,2.186269e+08,1.774119e+08,1.430176e+08,1.142823e+08,9.072941e+07,7.094795e+07,5.432124e+07,4.189605e+07,3.280686e+07,2.562972e+07,2.013453e+07,1.620424e+07,1.315460e+07,1.040263e+07,7.846039e+06,5.728553e+06,4.495501e+06,3.800058e+06,3.005911e+06,2.229857e+06,1.675695e+06,1.307405e+06,1.003372e+06,7.486895e+05,6.805768e+05,7.435860e+05,7.812756e+05,7.396855e+05,6.082531e+05,3.731120e+05,1.582579e+05,8.976509e+04,1.072953e+05,9.867716e+04,5.302092e+04,3.384872e+04,7.221184e+04,1.380859e+05,1.713877e+05,1.450725e+05,9.890296e+04,6.954686e+04,3.619885e+04,1.297526e+04,1.931888e+04,5.546763e+04,1.066047e+05,1.070156e+05,5.457797e+04,1.540227e+04,3.078500e+03,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00};       
   

  // flow rate in ml/min
  V_inj = 20. * TDatabase::ParamDB->UREA_INFLOW_SCALE / (3.0 * TDatabase::ParamDB->UREA_u_infty);
  // injection time
  V_inj *= TDatabase::ParamDB->UREA_inflow_time;
  V_inj = 1./V_inj;
  // scaling to units
  V_inj *= 60 * 1e6;

  Dscal(100,V_inj,inlet_f_L_seed);
        
  value = 0.;
      
    for (i=100;i>0;i--)
     { 
      if((L*L_max<=inlet_coord[i])&&(L*L_max>inlet_coord[i-1]))
      {
        m=(inlet_f_L_seed[i]-inlet_f_L_seed[i-1])/(inlet_coord[i] - inlet_coord[i-1]);
        value = m *(L*L_max- inlet_coord[i-1]);
        value += inlet_f_L_seed[i-1];
// 	if(L<1e-5)
//           OutPut(" L is !!!"<< L *L_max <<" i is !!!" <<i << " i is !!!"<<inlet_coord[i]<<endl);
        break;
       }
      }      
  
  return(value);
  
}// Inlet_PSD

// value of boundary condition
void BoundValue_psd(double x, double y, double z, double &value)
{
 double L = TDatabase::ParamDB->REACTOR_P29;
 
  value= 0.;   
  // inflow
  if( Is_InletBD(x, y, z) &&  (TDatabase::TimeDB->CURRENTTIME <= TDatabase::ParamDB->UREA_inflow_time) )
    {
     value = Inlet_PSD(L)/TDatabase::ParamDB->UREA_f_infty;
     //OutPut(L << " " << value << endl);      
    } 
    
}  

void BoundCondition_LminLMax(BoundCond &cond_Lmin, BoundCond &cond_Lmax)
{
  cond_Lmin = DIRICHLET;
//   cond_Lmin = NEUMANN;
  cond_Lmax = NEUMANN;
}

void BoundValue_LMin(double x, double y, double z, double *values)
 {
   
   cout<< "Nucleation implemented separately " << endl;
   exit(0);
   
   
   
   //nucleation
   // no nucleation, so f = 0 at L=0
   
//    double t;   
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
//   {
//     values[0] = 0.;
//     values[1] = 0.;
//     values[2] = 0.;
//   } 

//   double L_min = TDatabase::ParamDB->UREA_D_P_0;
//   double L_max = TDatabase::ParamDB->UREA_D_P_MAX;
//   double f_infty = TDatabase::ParamDB->UREA_f_infty;
//   double eps = 1e-8, val, a_min, g;
//   double f_in = TDatabase::ParamDB->UREA_f_infty;
// 



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
  int i;
  double eps, *coeff, *param;

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
      coeff[3] = param[2];  // u3      
//       cout<< "coeff[0] eps " << eps << " u1 " << param[0]  << "u2 " << param[1] << endl;
     }
    else
     {
      coeff[1] = 0;  // u1
      coeff[2] = 0;  // u2
      coeff[3] = 0;  // u2      
     }
     
    coeff[4] = 0;
    coeff[5] = 0.;
    coeff[6] = 0;  
  }
}


void BilinearCoeffs_Psd_Intl(int n_points, int N_Dim, double **Coords,
                             double **parameters, double **coeffs)
{
  int i;
  double eps, *coeff, b;
  double t = TDatabase::TimeDB->CURRENTTIME;

  b = -1e8;// negative, so that C will be taken from the PBS growth term
 
  if(TDatabase::ParamDB->REACTOR_P3)
    eps = 1.0/TDatabase::ParamDB->REACTOR_P3;
  else
    eps = 0.;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];  
   
    coeff[0] = eps;  // diffusion in L -direction    
    coeff[1] = b; // convection in L -direction
    coeff[2] = 0.;    // reaction L -direction
    coeff[3] = 0.;   // rhs in L -direction
  }
}

//input and output : dimless
double GetCsat(double C, double T)
{
 double val;
 double T_ref = TDatabase::ParamDB->REACTOR_P23;
 double C_ref = TDatabase::ParamDB->REACTOR_P25;
  
  val = (1.3045*(T_ref*T - 273.15) + 35.3642)/(C_ref*TDatabase::ParamDB->UREA_m_mol); 
  
  return(val);
}

double GetGrowth(double C, double T)
{
  double PbeGC = TDatabase::ParamDB->REACTOR_P27;
  double G, C_Sat, Mode = TDatabase::ParamDB->PBE_P0;

   C_Sat = GetCsat(C, T);

    // growth rate based on concentration
    if(C>C_Sat && (Mode==1 || Mode==4|| Mode==6) )
     { 
      G = PbeGC*sqrt((C - C_Sat )/C_Sat);       
      if(G<0.) G = 0.;
     }
    else
     { G = 0.; }

//           cout << "G: " << G << "C: " <<C <<"T: " << T <<endl;   
     
 return G;
 }
 
 
 // input and output in dimless form
void Get_GrowthAndB_Nuc(int N_Inputs, double *Inn, double *Out) 
{
//   double PbeGC = TDatabase::ParamDB->REACTOR_P27;
//   double C_Sat;
//   double T_ref = TDatabase::ParamDB->REACTOR_P23;
//   double C_ref = TDatabase::ParamDB->REACTOR_P25;
//   double alpha =  TDatabase::ParamDB->UREA_alfa_nuc;
//   double beta =  TDatabase::ParamDB->UREA_beta_nuc;
//   double f_ref = TDatabase::ParamDB->UREA_f_infty;
  double C, T, G, B_Nuc;
 
  if(N_Inputs!=2)
  {
    printf("N_Inputs != 2, needed Temp and conc to Get_GrowthAndB_Nuc: %d \n",N_Inputs); 
#ifdef _MPI
     MPI_Finalize();
#endif       
    exit(0);        
  }
  
  C = Inn[0];
  T = Inn[1];
  
//    C_Sat = (1.3045*(T_ref*T - 273.15) + 35.3642)/C_ref;

    // growth rate based on concentration
    G = GetGrowth(C, T); 
    // set constant
    B_Nuc = TDatabase::ParamDB->UREA_alfa_nuc;  
    if( G>0.)
     { B_Nuc /= (G*TDatabase::ParamDB->UREA_f_infty); }
    else
     { B_Nuc = 0.; } // no growth      
    
 
//     B_Nuc = C/C_Sat;
//   
//     if (B_Nuc != 0)
//      {
//       B_Nuc=log(B_Nuc);
//       B_Nuc *= B_Nuc;
//       B_Nuc = -beta/B_Nuc;
//       B_Nuc = alpha*exp(B_Nuc);  
//       
//       // set constant
// //       B_Nuc = TDatabase::ParamDB->UREA_alfa_nuc;    
// 
//       if( G>0.)
//        { B_Nuc /= (f_ref*G); }
//       else
//        { B_Nuc /= f_ref; } // no growth
// 
//        if(B_Nuc>1.)
//          B_Nuc = 1.;
// 
//      }
//     else
//      { B_Nuc = 0.0; } 
     
    Out[0] = G;
    Out[1] = B_Nuc;       
}

double PSD_BdValue(double x, double y, double z, double L, double c, double temp)
{
  double L_min = TDatabase::ParamDB->UREA_D_P_0;
  double L_max = TDatabase::ParamDB->UREA_D_P_MAX;
  double f_infty = TDatabase::ParamDB->UREA_f_infty;
  double eps = 1e-8, val, Inn[2], Out[2];
  double f_in = TDatabase::ParamDB->UREA_f_infty;

  // inflow
  if( Is_InletBD(x, y, z) )
   {  
    if(TDatabase::TimeDB->CURRENTTIME <= TDatabase::ParamDB->UREA_inflow_time)
     {       
      val = Inlet_PSD(L); // size of the seed crystals
      //OutPut(L << " " << val << endl);
      return(val);
     }
    else 
    { return(0.); }
   }

  L_min /= L_max;

   // at L_min
  if (fabs(L-L_min)<eps)
   {
    Inn[0] = c;
    Inn[1] = temp;    
    Get_GrowthAndB_Nuc(2, Inn, Out);
    val = Out[1];
    return(val);
   }
  return(0.0);
}

void GetExampleFileData(BoundCondFunct3D **BoundaryConditions, BoundValueFunct3D **BoundValues, 
                        DoubleFunct3D **InitiaValues, CoeffFct3D **BilinearCoeffs, 
                        int &N_PBEqns, int &N_IndepntScalarEqns, int *Disctypes)
{
  
   N_IndepntScalarEqns = 2;
   N_PBEqns = 1;
   #define __PBS__

   BilinearCoeffs[0] = BilinearCoeffs_Heat;
   BilinearCoeffs[1] = BilinearCoeffs_Conc;
   BilinearCoeffs[2] = BilinearCoeffs_Psd;

   
   BoundaryConditions[0] = BoundCondition_temp;
   BoundaryConditions[1] = BoundCondition_conc;
   BoundaryConditions[2] = BoundCondition_psd;

   
   BoundValues[0] = BoundValue_temp;
   BoundValues[1] = BoundValue_conc;
   BoundValues[2] = BoundValue_psd;

   
   InitiaValues[0] = InitialCondition_temp;
   InitiaValues[1] = InitialCondition_conc;
   InitiaValues[2] = InitialCondition_psd;

   
   Disctypes[0] = GALERKIN;
   Disctypes[1] = GALERKIN;
// //    Disctypes[2] = SDFEM;
   Disctypes[2] = GALERKIN;
}

void Generate1DMesh(TDomain *Domain, double Start, double End, int N_Cells)
{
  int i, j, N_V;
  int *Lines;
  double len, h, x, y, z;
  TVertex **Vetrex;
  TJoint *Joint;
  TBaseCell  **CellTree;

  
  N_V = N_Cells+1; 

//     double *X;
//   X = new double[N_V];
// 
//   for(i=0; i<N_V; i++)
//    X[i] = 1. + (1. - Start)*(tanh(2.75*(double(i)/double(N_Cells) - 1.)))/tanh(2.75);
// 
//   h = (End - Start)/N_Cells;
//   for(i=1; i<N_V; i++)
//    X[i] =  h*(double)i;
  
  
  if(N_Cells!=93)
   {
    printf("N_Cells_Internal != 112, Generate1DMesh : %d \n", N_Cells); 
#ifdef _MPI
     MPI_Finalize();
#endif       
    exit(0);
   } 

  
//   double X[113] = {0.001250000000000 , 0.004948299443979 , 0.006217671606001 , 0.007823152131089 , 
// 0.009849848192817 , 0.012405802335839 , 0.015627666201516 , 0.019687946143196 , 0.024804199643075 , 
// 0.028393316830700 , 0.031250666632101 , 0.035772828811153 , 0.039372952752786 , 0.042413137150944 , 
// 0.045070619540390 , 0.047446991851855 , 0.049606547407142 , 0.053436976294678 , 0.056785220390605 , 
// 0.059779281606086 , 0.062500166625532 , 0.065002773507971 , 0.067326227648337 , 0.069499512571455 , 
// 0.071544767303250 , 0.075317060474946 , 0.078745170560444 , 0.081898265580090 , 0.084825640943857 , 
// 0.087563813936174 , 0.090140678210433 , 0.092578078038834 , 0.094893477609017 , 0.097101084797404 , 
// 0.099212631825087 , 0.103185287607544 , 0.106873553596882 , 0.110323438897215 , 0.113570087453620 , 
// 0.116641021427289 , 0.119558244391400 , 0.122339657205393 , 0.125000041585273 , 0.130005277376165 , 
// 0.134652203946564 , 0.138998789266748 , 0.143089312023620 , 0.146958448032182 , 0.150633920105047 , 
// 0.154138281947656 , 0.157490157382677 , 0.163796361297538 , 0.169651123546876 , 0.175127479279515 , 
// 0.180281216202297 , 0.185156023145274 , 0.189786828693602 , 0.194202048758084 , 0.198425147902268 , 
// 0.206370468208258 , 0.213747007445228 , 0.220646784186760 , 0.227140086575032 , 0.233281959112378 , 
// 0.239116409077369 , 0.244679238288390 , 0.250000010253906 , 0.260010487342228 , 0.269304345055468 , 
// 0.277997519564343 , 0.286178568401423 , 0.293916843310072 , 0.301267789998808 , 0.308276515941197 , 
// 0.314980268830742 , 0.327592680129367 , 0.339302207508500 , 0.350254921410786 , 0.360562397349922 , 
// 0.370312013057423 , 0.379573625756073 , 0.388404067306966 , 0.396850266867541 , 0.412740909664792 , 
// 0.427493989953308 , 0.441293544971590 , 0.454280151067003 , 0.466563897289198 , 0.478232798228372 , 
// 0.489358457546173 , 0.500000002278646 , 0.520020957831925 , 0.538608674401517 , 0.555995024386395 , 
// 0.572357122891389 , 0.587833673431567 , 0.602535567444792 , 0.616553019893863 , 0.629960526177828 , 
// 0.655185349642304 , 0.678604405120686 , 0.700509833534509 , 0.721124785936175 , 0.740624017806563 , 
// 0.759147243604363 , 0.776808127061630 , 0.793700526500832 , 0.825481812641652 , 0.854987973672328 , 
// 0.882587084092698 , 0.908560296613240 , 0.956465591475152 , 1.000000000000000}; 

     //neu grid von ARAM, ref Bulk_3d4d.C
     double X[94] = {0, 0.000488281, 0.000615194, 0.000775097, 0.000976561, 0.00123039, 0.00155019, 0.00195312,0.00246078, 
                          0.00310039, 0.00390625, 0.00447154, 0.00492156, 0.0053016, 0.00563379, 0.00593084, 0.00620079,
                          0.00667959, 0.00709813, 0.00747238, 0.0078125, 0.00841576, 0.00894308, 0.00941462, 0.00984313, 
                          0.0106032, 0.0112676, 0.0118617, 0.0124016, 0.0133592, 0.0141962, 0.0149448, 0.015625, 0.0168315, 
                          0.0178862, 0.0188292, 0.0196863, 0.0204745, 0.0212064, 0.0218909, 0.0225352, 0.0231445, 0.0237233, 
                          0.0242752, 0.0248031, 0.0267184, 0.0283925, 0.0298896, 0.03125, 0.0336631, 0.0357723, 0.0393725, 
                          0.0450703, 0.0496062, 0.0534368, 0.056785, 0.0597791, 0.0625, 0.0673261, 0.0715446, 0.0753169,
                          0.0787451, 0.0818982, 0.0848255, 0.0875637, 0.0901406, 0.0948934, 0.0992125, 0.106873, 0.11357,
                          0.119558, 0.125, 0.134652, 0.143089, 0.150634, 0.15749, 0.169651, 0.180281, 0.189787, 0.198425, 
                          0.213747, 0.22714, 0.239117, 0.25, 0.269304, 0.286179, 0.301268, 0.31498, 0.360562, 0.39685, 
                          0.5, 0.629961, 0.793701, 1};


//   X[0] = Start;
//   X[N_V-1] = End;

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



